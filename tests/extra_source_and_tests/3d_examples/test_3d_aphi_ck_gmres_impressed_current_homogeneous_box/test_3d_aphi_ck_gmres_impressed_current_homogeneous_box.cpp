/**
 * Stage 9C-1: homogeneous box with impressed current source (physical RHS, not MMS).
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
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 10.0;
    const Real tolerance = 1.0e-5;
    const Real impressed_current_amplitude = 5.0;
    const Real source_exterior_ratio_min = 1.2;
    const Real min_source_region_solution_block_max = 1.0e-3;

    AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};
    const Vecd current_real(0.0, 0.0, 1.0);
    const Vecd current_imag(0.0, 0.0, 0.1);

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignImpressedCurrentRhsCK> assign_impressed_current(
        test_body.body, names.rhs, source_region, current_real, current_imag, impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    AphiMatrixFreeSolverOptions solver_options = defaultMatrixFreeGMRESOptions(tolerance);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    initialize_aphi_variables.exec();
    set_material.exec();
    assign_impressed_current.exec();
    test_body.updateRelations();

    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult result = solver.solve();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real source_region_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, body_length,
                                    body_height, body_width, core_shell, source_region);
    const Real exterior_solution_block_max =
        coreBlockMaxAbsNormOutsideRegion(particles, names.solution, positions, total_real_particles, body_length,
                                         body_height, body_width, core_shell, source_region);
    const Real source_exterior_ratio = source_region_solution_block_max / (exterior_solution_block_max + TinyReal);

    const bool passed = gmresConvergencePassed(result, tolerance) &&
                        source_region_solution_block_max >= min_source_region_solution_block_max &&
                        source_exterior_ratio >= source_exterior_ratio_min;

    std::cout << "test_3d_aphi_ck_gmres_impressed_current_homogeneous_box"
              << " particles=" << total_real_particles << " init_norm=" << result.initial_residual_norm
              << " outer=" << result.outer_iteration_count << " arnoldi=" << result.arnoldi_step_count
              << " rel_res=" << result.final_relative_residual << " true_rel_res=" << result.final_true_relative_residual
              << " monotonic_outer=" << (result.monotonic_outer_residual ? 1 : 0)
              << " source_region_solution_block_max=" << source_region_solution_block_max
              << " exterior_solution_block_max=" << exterior_solution_block_max
              << " source_exterior_ratio=" << source_exterior_ratio
              << " breakdown_code=" << result.breakdown_code_name << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
