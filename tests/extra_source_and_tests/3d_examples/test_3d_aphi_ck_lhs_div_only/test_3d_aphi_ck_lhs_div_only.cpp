#include "electromagnetic_dynamics/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignDivOnlyFieldsCK : public LocalDynamics
{
  public:
    explicit AssignDivOnlyFieldsCK(SPHBody &sph_body, const AphiVariableNames &names)
        : LocalDynamics(sph_body),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(names.solution.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(names.solution.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(names.solution.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(names.solution.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            a_real_[index_i] = position;
            a_imag_[index_i] = Vecd(position[1], position[2], position[0]);
            phi_real_[index_i] = 0.0;
            phi_imag_[index_i] = 0.0;
        }

      protected:
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real sigma = 1.0;
    const Real omega = 1.0;
    const Real validation_threshold = 0.35;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width);
    test_body.sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, 1.0, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, 1.0, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivOnlyFieldsCK> assign_fields(test_body.body, names);

    AphiLhsAssemblyOptions options;
    options.terms.laplace_a = false;
    options.terms.laplace_phi = false;
    options.terms.reaction = false;
    options.terms.grad_phi_coupling = false;
    options.terms.div_sigma_a_coupling = true;
    options.div_mode = AphiDivCouplingMode::PairwiseUncorrected;
    options.omega = omega;

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble(test_body.body, test_body.inner_ck, names, options);

    initialize_aphi_variables.exec();
    set_material.exec();
    assign_fields.exec();
    test_body.updateRelations();
    assemble.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *lhs_phi_real = particles.getVariableDataByName<Real>(names.lhs.phi_real);
    const Real *lhs_phi_imag = particles.getVariableDataByName<Real>(names.lhs.phi_imag);

    const Real expected_phi_real = 0.0;
    const Real expected_phi_imag = -omega * sigma * 3.0;

    Real core_max_error = 0.0;
    size_t core_particles = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        ++core_particles;
        core_max_error = std::max(core_max_error, std::abs(lhs_phi_real[i] - expected_phi_real));
        core_max_error = std::max(core_max_error, std::abs(lhs_phi_imag[i] - expected_phi_imag));
    }

    const bool passed = core_particles > 0 && core_max_error < validation_threshold;

    std::cout << "test_3d_aphi_ck_lhs_div_only"
              << " core_particles=" << core_particles
              << " core_max_error=" << core_max_error
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
