#include "electromagnetic_dynamics/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignPolynomialAphiFieldsCK : public LocalDynamics
{
  public:
    explicit AssignPolynomialAphiFieldsCK(SPHBody &sph_body, const AphiVariableNames &names)
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
            const Real x = position[0];
            const Real y = position[1];
            const Real z = position[2];
            const Real quadratic_sum = x * x + y * y + z * z;

            phi_real_[index_i] = quadratic_sum;
            phi_imag_[index_i] = 0.5 * quadratic_sum;
            a_real_[index_i] = Vecd(x * x, y * y, z * z);
            a_imag_[index_i] = Vecd(y * y, z * z, x * x);
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

Real coreBlockMaxAbsLhs(BaseParticles &particles, const AphiVariableNames &names, const Vecd *positions,
                        size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                        Real core_shell)
{
    const Vecd *lhs_a_real = particles.getVariableDataByName<Vecd>(names.lhs.a_real);
    const Vecd *lhs_a_imag = particles.getVariableDataByName<Vecd>(names.lhs.a_imag);
    const Real *lhs_phi_real = particles.getVariableDataByName<Real>(names.lhs.phi_real);
    const Real *lhs_phi_imag = particles.getVariableDataByName<Real>(names.lhs.phi_imag);
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        max_value = std::max(max_value, lhs_a_real[i].norm());
        max_value = std::max(max_value, lhs_a_imag[i].norm());
        max_value = std::max(max_value, std::abs(lhs_phi_real[i]));
        max_value = std::max(max_value, std::abs(lhs_phi_imag[i]));
    }
    return max_value;
}

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
    const Real nu = 1.0;
    const Real omega = 0.5;
    const Real min_lhs_magnitude = 1.0e-3;
    const Real max_lhs_magnitude = 1.0e6;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width);
    test_body.sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignPolynomialAphiFieldsCK> assign_fields(test_body.body, names);

    AphiLhsAssemblyOptions options;
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

    size_t core_particles = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            ++core_particles;
        }
    }

    const Real core_max_lhs = coreBlockMaxAbsLhs(particles, names, positions, total_real_particles, body_length,
                                               body_height, body_width, core_shell);
    const bool passed = core_particles > 0 && core_max_lhs > min_lhs_magnitude && core_max_lhs < max_lhs_magnitude;

    std::cout << "test_3d_aphi_ck_lhs_polynomial_manufactured"
              << " core_particles=" << core_particles
              << " core_max_lhs=" << core_max_lhs
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
