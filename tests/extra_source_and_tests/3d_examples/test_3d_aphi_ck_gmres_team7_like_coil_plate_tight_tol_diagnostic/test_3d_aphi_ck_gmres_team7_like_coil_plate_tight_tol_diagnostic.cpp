/**
 * Stage 9C-3 diagnostic: TEAM7-like case at target tol=1e-5 (informational PC floor record).
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
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const Real target_tolerance = 1.0e-5;
    const UnsignedInt restart_dimension = 30;
    const UnsignedInt max_outer_iterations = 100;

    const AphiTeam7LikeUnitBoxLayout layout;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, AssignTeam7LikeRegionMaterialsCK> assign_material(test_body.body, layout,
                                                                                         names.material);
    StateDynamics<MainExecutionPolicy, AssignImpressedCurrentRhsCK> assign_coil_source(
        test_body.body, names.rhs, layout.coil, coil_current_real, coil_current_imag, impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(target_tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();

    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult result = solver.solve();
    const bool converged_at_target = gmresConvergencePassed(result, target_tolerance);

    std::cout << "test_3d_aphi_ck_gmres_team7_like_coil_plate_tight_tol_diagnostic"
              << " target_tolerance=" << target_tolerance << " init_norm=" << result.initial_residual_norm
              << " outer=" << result.outer_iteration_count << " arnoldi=" << result.arnoldi_step_count
              << " rel_res=" << result.final_relative_residual << " true_rel_res=" << result.final_true_relative_residual
              << " converged_at_1e-5=" << (converged_at_target ? 1 : 0)
              << " breakdown_code=" << result.breakdown_code_name << " passed=1" << std::endl;

    return 0;
}
