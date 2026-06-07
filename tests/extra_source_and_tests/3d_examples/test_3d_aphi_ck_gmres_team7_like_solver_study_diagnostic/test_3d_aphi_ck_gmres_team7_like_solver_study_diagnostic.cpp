/**
 * Stage 9D-0 diagnostic: TEAM7-like full-contrast GMRES study.
 * - restart sweep m=30/50/80 with block-Jacobi PC
 * - no-PC vs block-Jacobi at m=30
 * - component-wise and region-wise true residual breakdown
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <array>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

void logGMRESResult(std::ostream &out, const char *label, const AphiGMRESResult &result)
{
    out << " " << label << " outer=" << result.outer_iteration_count << " arnoldi=" << result.arnoldi_step_count
        << " rel_res=" << result.final_relative_residual << " true_rel_res=" << result.final_true_relative_residual
        << " converged=" << (result.converged ? 1 : 0)
        << " breakdown_code=" << AphiGMRESBreakdownCodeName(result.breakdown_code);
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 3.0 * dp_0;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const Real tolerance = 1.0e-5;
    const UnsignedInt max_workspace_dimension = 80;
    const UnsignedInt max_outer_iterations = 100;
    const std::array<UnsignedInt, 3> restart_values = {30, 50, 80};

    const AphiTeam7LikeUnitBoxLayout layout;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    const AphiGMRESWorkspaceNames workspace = buildAphiGMRESWorkspaceNames(max_workspace_dimension);
    RegisterAphiGMRESWorkspaceCK register_gmres_workspace(test_body.body, max_workspace_dimension);

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

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    std::cout << "test_3d_aphi_ck_gmres_team7_like_solver_study_diagnostic tolerance=" << tolerance
              << " max_outer=" << max_outer_iterations;

    AphiGMRESResult breakdown_reference_result;
    Real breakdown_reference_norm = 0.0;

    std::cout << " restart_sweep_pc=1";
    for (const UnsignedInt restart_dimension : restart_values)
    {
        const AphiGMRESResult result = runTeam7LikeGMRESDirect(
            test_body, names, workspace, options, zero_solution, tolerance, restart_dimension, max_outer_iterations,
            true);
        logGMRESResult(std::cout, ("m" + std::to_string(restart_dimension) + "_bj").c_str(), result);
    }

    std::cout << " pc_comparison_m30";
    const AphiGMRESResult no_pc_result = runTeam7LikeGMRESDirect(
        test_body, names, workspace, options, zero_solution, tolerance, 30, max_outer_iterations, false,
        AphiBlockJacobiPreconditionerKind::Decoupled);
    logGMRESResult(std::cout, "no_pc", no_pc_result);

    const AphiGMRESResult decoupled_result = runTeam7LikeGMRESDirect(
        test_body, names, workspace, options, zero_solution, tolerance, 30, max_outer_iterations, true,
        AphiBlockJacobiPreconditionerKind::Decoupled);
    logGMRESResult(std::cout, "decoupled_m30", decoupled_result);

    const AphiGMRESResult block_jacobi_result = runTeam7LikeGMRESDirect(
        test_body, names, workspace, options, zero_solution, tolerance, 30, max_outer_iterations, true,
        AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8);
    logGMRESResult(std::cout, "coupled_m30", block_jacobi_result);

    logTrueResidualDiagnostics(std::cout, particles, names, positions, total_real_particles, body_length, body_height,
                               body_width, core_shell, block_jacobi_result.initial_residual_norm, layout,
                               "residual_breakdown");

    std::cout << " passed=1" << std::endl;
    return 0;
}
