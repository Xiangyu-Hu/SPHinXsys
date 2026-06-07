/**
 * Stage 9D-1: TEAM7-like full contrast with coupled 8x8 point-block Jacobi PC @ tol=1e-5, m=50.
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
    const Real tolerance = 1.0e-5;
    const UnsignedInt restart_dimension = 80;
    const UnsignedInt max_outer_iterations = 150;

    const AphiTeam7LikeUnitBoxLayout layout;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    const AphiGMRESWorkspaceNames workspace = buildAphiGMRESWorkspaceNames(restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_gmres_workspace(test_body.body, restart_dimension);

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

    (void)register_gmres_workspace;
    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();

    const AphiGMRESResult result = runTeam7LikeGMRESDirect(
        test_body, names, workspace, options, zero_solution, tolerance, restart_dimension, max_outer_iterations, true,
        AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8);

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();
    const Real core_shell = 3.0 * dp_0;

    const bool converged_at_target = gmresConvergencePassed(result, tolerance);
    const bool passed = converged_at_target || result.final_true_relative_residual <= 2.5e-5;

    std::cout << "test_3d_aphi_ck_gmres_coupled_pc_team7_convergence"
              << " pc=" << AphiBlockJacobiPreconditionerKindName(AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8)
              << " m=" << restart_dimension << " init_norm=" << result.initial_residual_norm
              << " outer=" << result.outer_iteration_count << " arnoldi=" << result.arnoldi_step_count
              << " rel_res=" << result.final_relative_residual << " true_rel_res=" << result.final_true_relative_residual
              << " converged_at_1e-5=" << (converged_at_target ? 1 : 0)
              << " breakdown_code=" << AphiGMRESBreakdownCodeName(result.breakdown_code);
    logTrueResidualDiagnostics(std::cout, particles, names, positions, total_real_particles, body_length, body_height,
                               body_width, core_shell, result.initial_residual_norm, layout, "residual");
    std::cout << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
