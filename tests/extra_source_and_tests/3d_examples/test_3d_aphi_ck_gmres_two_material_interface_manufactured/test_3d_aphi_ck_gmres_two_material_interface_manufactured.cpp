/**
 * Stage 9C-2: two-material conductor/air interface with discrete MMS (variable sigma).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real sigma_conductor = 2.0;
    const Real sigma_air = 1.0e-4;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 10.0;
    const Real tolerance = 1.0e-5;
    const Real max_discrete_mms_defect_factor = 10.0;
    const Real x_interface = 0.5;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, sigma_conductor, nu, names);
    StateDynamics<MainExecutionPolicy, AssignPiecewiseSigmaHalfSpaceCK> assign_material(
        test_body.body, x_interface, sigma_conductor, sigma_air, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiManufacturedFieldsCK> assign_exact(test_body.body,
                                                                                             names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact_reference(test_body.body, names.r_hat, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_exact(test_body.body, test_body.inner(), names.solution, names.lhs,
                                                             names.material, omega, options);
    AphiMatrixFreeSolverOptions solver_options = defaultMatrixFreeGMRESOptions(tolerance);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_exact.exec();
    test_body.updateRelations();

    apply_exact.exec();
    copy_exact_reference.exec();
    copy_lhs_to_rhs.exec();
    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult result = solver.solve();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real continuous_field_relative_error =
        coreRelativeBlockError(particles, names.solution, names.r_hat, positions, total_real_particles, body_length,
                               body_height, body_width, core_shell);
    const Real discrete_mms_defect =
        coreDiscreteMmsOperatorDefect(test_body, names, options, names.solution, names.r_hat);

    const bool passed = gmresConvergencePassed(result, tolerance) &&
                        discrete_mms_defect <= max_discrete_mms_defect_factor * tolerance;

    std::cout << "test_3d_aphi_ck_gmres_two_material_interface_manufactured"
              << " particles=" << total_real_particles << " init_norm=" << result.initial_residual_norm
              << " outer=" << result.outer_iteration_count << " arnoldi=" << result.arnoldi_step_count
              << " rel_res=" << result.final_relative_residual << " true_rel_res=" << result.final_true_relative_residual
              << " discrete_mms_defect=" << discrete_mms_defect
              << " continuous_field_relative_error=" << continuous_field_relative_error
              << " monotonic_outer=" << (result.monotonic_outer_residual ? 1 : 0)
              << " breakdown_code=" << result.breakdown_code_name << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
