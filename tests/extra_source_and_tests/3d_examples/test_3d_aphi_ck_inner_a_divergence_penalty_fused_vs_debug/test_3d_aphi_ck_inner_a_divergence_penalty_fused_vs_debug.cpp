/** Stage 10.7 inner: fused apply vs debug LHS with A-divergence penalty enabled. */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignSeparableAphiFieldsCK : public LocalDynamics
{
  public:
    explicit AssignSeparableAphiFieldsCK(SPHBody &sph_body, const AphiBlockNames &block_names)
        : LocalDynamics(sph_body),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag))
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
            const Real x = position[0];
            const Real y = position[1];
            const Real z = position[2];
            const Real pi = Pi;

            a_real_[index_i] = Vecd(std::sin(pi * x), std::sin(pi * y), std::sin(pi * z));
            a_imag_[index_i] = Vecd(std::cos(pi * x), std::cos(pi * y), std::cos(pi * z));
            phi_real_[index_i] = std::sin(pi * x) * std::cos(pi * y);
            phi_imag_[index_i] = std::cos(pi * x) * std::sin(pi * z);
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
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real a_divergence_penalty = 50.0;
    const Real cross_path_threshold = 5.0e-3;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_fields(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = true;
    options.a_divergence_penalty = a_divergence_penalty;

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble_debug(test_body.body, test_body.inner(), names, options);
    AphiApplyDynamicsBundle<MainExecutionPolicy> fused_apply(test_body.body, test_body.inner(), names.solution, names.v,
                                                             names.material, omega, options);

    initialize_aphi_variables.exec();
    set_material.exec();
    assign_fields.exec();
    test_body.updateRelations();

    assemble_debug.exec();
    fused_apply.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    const Real core_max_difference =
        coreBlockMaxAbsDifference(particles, names.lhs, names.v, positions, total_real_particles, body_length,
                                  body_height, body_width, core_shell);
    const bool passed = core_max_difference < cross_path_threshold;

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_fused_vs_debug"
              << " core_max_difference=" << core_max_difference
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
