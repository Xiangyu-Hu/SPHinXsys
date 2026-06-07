#ifndef APHI_GMRES_ROBUSTNESS_SWEEP_HELPERS_H
#define APHI_GMRES_ROBUSTNESS_SWEEP_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiManufacturedRobustnessParams
{
    Real dp_0 = 0.1;
    Real sigma = 2.0;
    Real nu = 1.5;
    Real omega = 1.25;
    Real phi_gauge_penalty = 100.0;
};

struct AphiManufacturedRobustnessResult
{
    size_t total_real_particles = 0;
    Real initial_residual_norm = 0.0;
    Real final_relative_residual = 0.0;
    Real final_true_relative_residual = 0.0;
    UnsignedInt outer_iteration_count = 0;
    UnsignedInt arnoldi_step_count = 0;
    bool monotonic_outer_residual = true;
    bool converged = false;
    bool breakdown = false;
    const char *breakdown_code_name = "None";
};

class AssignSeparableAphiManufacturedFieldsCK : public LocalDynamics
{
  public:
    explicit AssignSeparableAphiManufacturedFieldsCK(SPHBody &sph_body, const AphiBlockNames &block_names)
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

inline AphiManufacturedRobustnessResult runFullAphiManufacturedGMRES(int ac, char *av[],
                                                                     const AphiManufacturedRobustnessParams &params,
                                                                     Real tolerance, UnsignedInt restart_dimension,
                                                                     UnsignedInt max_outer_iterations)
{
    AphiManufacturedRobustnessResult sweep_result;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * params.dp_0;

    AphiLhsTestBody test_body(params.dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, params.sigma,
                                                                                        params.nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, params.sigma, params.nu,
                                                                               names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiManufacturedFieldsCK> assign_exact(test_body.body,
                                                                                             names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = params.omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = params.phi_gauge_penalty;

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_exact(test_body.body, test_body.inner(), names.solution, names.lhs,
                                                             names.material, params.omega, options);
    AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    initialize_aphi_variables.exec();
    set_material.exec();
    assign_exact.exec();
    test_body.updateRelations();

    apply_exact.exec();
    copy_lhs_to_rhs.exec();
    zero_solution.exec();
    test_body.updateRelations();

    sweep_result.total_real_particles = test_body.body.getBaseParticles().TotalRealParticles();

    const AphiMatrixFreeSolverResult result = solver.solve();
    sweep_result.initial_residual_norm = result.initial_residual_norm;
    sweep_result.final_relative_residual = result.final_relative_residual;
    sweep_result.final_true_relative_residual = result.final_true_relative_residual;
    sweep_result.outer_iteration_count = result.outer_iteration_count;
    sweep_result.arnoldi_step_count = result.arnoldi_step_count;
    sweep_result.monotonic_outer_residual = result.monotonic_outer_residual;
    sweep_result.converged = result.converged;
    sweep_result.breakdown = result.breakdown;
    sweep_result.breakdown_code_name = result.breakdown_code_name;
    return sweep_result;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_ROBUSTNESS_SWEEP_HELPERS_H
