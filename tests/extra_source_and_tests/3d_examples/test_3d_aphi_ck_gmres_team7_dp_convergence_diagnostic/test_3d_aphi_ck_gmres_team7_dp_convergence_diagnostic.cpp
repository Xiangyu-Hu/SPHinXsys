/**
 * Stage 9E-3: TEAM7-like unit-box dp convergence sweep (informational diagnostic).
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

struct Team7DpSweepResult
{
    Real dp = 0.0;
    size_t particles = 0;
    Real init_norm = 0.0;
    UnsignedInt outer = 0;
    UnsignedInt arnoldi = 0;
    Real rel_res = 0.0;
    Real true_rel_res = 0.0;
    Real conductor_solution_block_max = 0.0;
    Real coil_solution_block_max = 0.0;
    bool converged = false;
};

Team7DpSweepResult runTeam7DpCase(int ac, char *av[], Real dp_0)
{
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const Real tolerance = 5.0e-4;
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;

    const AphiTeam7LikeUnitBoxLayout layout;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);

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
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();

    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult result = solver.solve();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    Team7DpSweepResult sweep;
    sweep.dp = dp_0;
    sweep.particles = total_real_particles;
    sweep.init_norm = result.initial_residual_norm;
    sweep.outer = result.outer_iteration_count;
    sweep.arnoldi = result.arnoldi_step_count;
    sweep.rel_res = result.final_relative_residual;
    sweep.true_rel_res = result.final_true_relative_residual;
    sweep.conductor_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, body_length,
                                    body_height, body_width, core_shell, layout.conductor);
    sweep.coil_solution_block_max = coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles,
                                                       body_length, body_height, body_width, core_shell, layout.coil);
    sweep.converged = gmresConvergencePassed(result, tolerance);
    return sweep;
}

void printSweepResult(const Team7DpSweepResult &result)
{
    std::cout << " dp=" << result.dp << " particles=" << result.particles << " init_norm=" << result.init_norm
              << " outer=" << result.outer << " arnoldi=" << result.arnoldi << " rel_res=" << result.rel_res
              << " true_rel_res=" << result.true_rel_res << " conductor_solution_block_max=" << result.conductor_solution_block_max
              << " coil_solution_block_max=" << result.coil_solution_block_max << " converged=" << (result.converged ? 1 : 0);
}

} // namespace

int main(int ac, char *av[])
{
    const std::array<Real, 3> dp_values = {0.2, 0.1, 0.075};

    std::cout << "test_3d_aphi_ck_gmres_team7_dp_convergence_diagnostic";
    for (const Real dp : dp_values)
    {
        const Team7DpSweepResult result = runTeam7DpCase(ac, av, dp);
        printSweepResult(result);
    }
    std::cout << " passed=1" << std::endl;
    return 0;
}
