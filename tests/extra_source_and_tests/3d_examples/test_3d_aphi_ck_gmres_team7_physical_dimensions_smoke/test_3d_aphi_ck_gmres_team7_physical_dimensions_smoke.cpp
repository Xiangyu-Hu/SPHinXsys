/**
 * Stage 9E-2: TEAM7 canonical physical box (1.2 x 1.0 x 0.3 m) with fraction-mapped layout.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <algorithm>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = AphiTeam7PhysicalDimensions::length;
    const Real body_height = AphiTeam7PhysicalDimensions::height;
    const Real body_width = AphiTeam7PhysicalDimensions::width;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell =
        std::min(Real(2.0) * dp_0, Real(0.125) * std::min({body_length, body_height, body_width}));
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const Real tolerance = 5.0e-4;
    const Real min_conductor_solution_block_max = 1.0e-4;
    const Real min_coil_solution_block_max = 1.0e-4;
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;

    const AphiTeam7LikeUnitBoxLayout layout =
        buildTeam7LayoutForBox(body_length, body_height, body_width);
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
    const Real conductor_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, body_length,
                                    body_height, body_width, core_shell, layout.conductor);
    const Real coil_solution_block_max = coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles,
                                                            body_length, body_height, body_width, core_shell,
                                                            layout.coil);

    const bool passed = gmresConvergencePassed(result, tolerance) &&
                        conductor_solution_block_max >= min_conductor_solution_block_max &&
                        coil_solution_block_max >= min_coil_solution_block_max;

    std::cout << "test_3d_aphi_ck_gmres_team7_physical_dimensions_smoke"
              << " body=" << body_length << "x" << body_height << "x" << body_width << " dp=" << dp_0
              << " particles=" << total_real_particles << " init_norm=" << result.initial_residual_norm
              << " outer=" << result.outer_iteration_count << " arnoldi=" << result.arnoldi_step_count
              << " rel_res=" << result.final_relative_residual << " true_rel_res=" << result.final_true_relative_residual
              << " conductor_solution_block_max=" << conductor_solution_block_max
              << " coil_solution_block_max=" << coil_solution_block_max
              << " monotonic_outer=" << (result.monotonic_outer_residual ? 1 : 0)
              << " breakdown_code=" << result.breakdown_code_name << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
