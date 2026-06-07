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
    explicit AssignSeparableAphiFieldsCK(SPHBody &sph_body, const AphiVariableNames &names)
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
    const Real validation_threshold = 1.0e-8;
    /** Cross-path (debug multi-kernel vs fused single-kernel) SYCL FP tolerance. */
    const Real cross_path_threshold = 1.0e-4;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_fields(test_body.body, names);

    AphiLhsAssemblyOptions options;
    options.omega = omega;

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble_debug(test_body.body, test_body.inner(), names, options);
    AphiApplyDynamicsBundle<MainExecutionPolicy> fused_apply(test_body.body, test_body.inner(), names.solution, names.v,
                                                             names.material, omega, options);
    AphiApplyDynamicsBundle<MainExecutionPolicy> fused_apply_repeat(test_body.body, test_body.inner(), names.solution,
                                                                    names.s, names.material, omega, options);

    initialize_aphi_variables.exec();
    set_material.exec();
    assign_fields.exec();
    test_body.updateRelations();

    assemble_debug.exec();
    fused_apply.exec();
    fused_apply_repeat.exec();

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

    const Real core_max_difference =
        coreBlockMaxAbsDifference(particles, names.lhs, names.v, positions, total_real_particles, body_length,
                                  body_height, body_width, core_shell);
    const Real core_fused_repeat_difference =
        coreBlockMaxAbsDifference(particles, names.v, names.s, positions, total_real_particles, body_length,
                                  body_height, body_width, core_shell);
    const bool passed = core_particles > 0 && core_max_difference < cross_path_threshold &&
                        core_fused_repeat_difference < validation_threshold;

    std::cout << "test_3d_aphi_ck_fused_apply_vs_debug_lhs"
              << " core_particles=" << core_particles
              << " core_max_difference=" << core_max_difference
              << " core_fused_repeat_difference=" << core_fused_repeat_difference
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
