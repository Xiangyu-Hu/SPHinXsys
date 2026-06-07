/**
 * Stage 7: constant-field nullspace diagnostic for pure Laplace.
 * Reports ||K(constant)|| to confirm Laplace-only systems are singular near constants.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignConstantFieldsCK : public LocalDynamics
{
  public:
    explicit AssignConstantFieldsCK(SPHBody &sph_body, const AphiBlockNames &block_names, Real phi_value,
                                  const Vecd &a_value)
        : LocalDynamics(sph_body), phi_value_(phi_value), a_value_(a_value),
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
            : phi_value_(encloser.phi_value_), a_value_(encloser.a_value_),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            a_real_[index_i] = a_value_;
            a_imag_[index_i] = a_value_;
            phi_real_[index_i] = phi_value_;
            phi_imag_[index_i] = phi_value_;
        }

      protected:
        Real phi_value_;
        Vecd a_value_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    Real phi_value_;
    Vecd a_value_;
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
    const Real constant_phi = 1.0;
    const Vecd constant_a(1.0, 1.0, 1.0);
    const Real validation_threshold = 1.0e-5;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantFieldsCK> assign_constant(test_body.body, names.solution, constant_phi,
                                                                               constant_a);

    AphiLhsAssemblyOptions phi_options;
    phi_options.terms.laplace_phi = true;
    phi_options.terms.laplace_a = false;
    phi_options.terms.reaction = false;
    phi_options.terms.grad_phi_coupling = false;
    phi_options.terms.div_sigma_a_coupling = false;

    AphiLhsAssemblyOptions a_options = phi_options;
    a_options.terms.laplace_phi = false;
    a_options.terms.laplace_a = true;

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_phi(test_body.body, test_body.inner(), names.solution, names.lhs,
                                                           names.material, phi_options.omega, phi_options);
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_a(test_body.body, test_body.inner(), names.solution, names.lhs,
                                                         names.material, a_options.omega, a_options);

    initialize_aphi_variables.exec();
    set_material.exec();
    assign_constant.exec();
    test_body.updateRelations();

    apply_phi.exec();
    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real phi_constant_norm = coreBlockOperatorNorm(particles, names.lhs, positions, total_real_particles, body_length,
                                                         body_height, body_width, core_shell);

    apply_a.exec();
    const Real a_constant_norm = coreBlockOperatorNorm(particles, names.lhs, positions, total_real_particles, body_length,
                                                       body_height, body_width, core_shell);

    /** K(const)=0 confirms constant lies in Laplace nullspace (singular operator). */
    const bool passed = phi_constant_norm < validation_threshold && a_constant_norm < validation_threshold;

    std::cout << "test_3d_aphi_ck_laplace_constant_nullspace_diagnostic"
              << " core_K_constant_phi_norm=" << phi_constant_norm
              << " core_K_constant_a_norm=" << a_constant_norm
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
