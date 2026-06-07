#ifndef APHI_INNER_A_DIVERGENCE_PENALTY_SWEEP_HELPERS_H
#define APHI_INNER_A_DIVERGENCE_PENALTY_SWEEP_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiInnerADivergencePenaltySweepRow
{
    Real a_divergence_penalty = 0.0;
    bool use_a_divergence_penalty = false;
    bool converged = false;
    UnsignedInt outer_iterations = 0;
    UnsignedInt restart_dimension = 0;
    UnsignedInt max_outer_iterations = 0;
    Real final_relative_residual = 0.0;
    Real final_true_relative_residual = 0.0;
    Real final_res_a_fraction = 0.0;
    Real final_res_phi_fraction = 0.0;
    Real source_region_solution_block_max = 0.0;
    Real source_exterior_ratio = 0.0;
    AphiDivAReductionMetrics global_div_a{};
    AphiDivAReductionMetrics source_region_div_a{};
};

inline AphiInnerADivergencePenaltySweepRow runInnerADivergencePenaltySweepRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real a_divergence_penalty,
    bool use_a_divergence_penalty, Real tolerance, const benchmark::AphiBoxRegion &source_region,
    const Vecd &current_real, const Vecd &current_imag, Real impressed_current_amplitude,
    UnsignedInt restart_dimension = 50, UnsignedInt max_outer_iterations = 100)
{
    AphiInnerADivergencePenaltySweepRow row;
    row.a_divergence_penalty = a_divergence_penalty;
    row.use_a_divergence_penalty = use_a_divergence_penalty;
    row.restart_dimension = restart_dimension;
    row.max_outer_iterations = max_outer_iterations;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignImpressedCurrentRhsCK> assign_impressed_current(
        test_body.body, names.rhs, source_region, current_real, current_imag, impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = use_a_divergence_penalty;
    options.a_divergence_penalty = a_divergence_penalty;

    const AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    initialize_aphi.exec();
    set_material.exec();
    assign_impressed_current.exec();
    test_body.updateRelations();
    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult result = solver.solve();
    row.converged = result.converged;
    row.outer_iterations = result.outer_iteration_count;
    row.final_relative_residual = result.final_relative_residual;
    row.final_true_relative_residual = result.final_true_relative_residual;
    row.final_res_a_fraction = result.final_res_a_fraction;
    row.final_res_phi_fraction = result.final_res_phi_fraction;

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    row.source_region_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, body_length,
                                     body_height, body_width, core_shell, source_region);
    const Real exterior_solution_block_max =
        coreBlockMaxAbsNormOutsideRegion(particles, names.solution, positions, total_real_particles, body_length,
                                         body_height, body_width, core_shell, source_region);
    row.source_exterior_ratio = row.source_region_solution_block_max / (exterior_solution_block_max + TinyReal);

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names);
    row.global_div_a = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles);
    row.source_region_div_a = hostBodyDivAReductionMetrics(
        particles, names, positions, total_real_particles,
        [&](const Vecd &position) { return benchmark::insideBoxRegion(position, source_region); });
    return row;
}

inline AphiInnerADivergencePenaltySweepRow runInnerSolenoidalADivergencePenaltySweepRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real a_divergence_penalty,
    bool use_a_divergence_penalty, Real tolerance, const benchmark::AphiBoxRegion &source_region,
    Real solenoidal_current_amplitude, Real imag_to_real_ratio = 0.1, UnsignedInt restart_dimension = 50,
    UnsignedInt max_outer_iterations = 100)
{
    AphiInnerADivergencePenaltySweepRow row;
    row.a_divergence_penalty = a_divergence_penalty;
    row.use_a_divergence_penalty = use_a_divergence_penalty;
    row.restart_dimension = restart_dimension;
    row.max_outer_iterations = max_outer_iterations;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignSolenoidalCurlCurrentRhsCK> assign_solenoidal_current(
        test_body.body, names.rhs, source_region, solenoidal_current_amplitude, imag_to_real_ratio);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = use_a_divergence_penalty;
    options.a_divergence_penalty = a_divergence_penalty;

    const AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    initialize_aphi.exec();
    set_material.exec();
    assign_solenoidal_current.exec();
    test_body.updateRelations();
    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult result = solver.solve();
    row.converged = result.converged;
    row.outer_iterations = result.outer_iteration_count;
    row.final_relative_residual = result.final_relative_residual;
    row.final_true_relative_residual = result.final_true_relative_residual;
    row.final_res_a_fraction = result.final_res_a_fraction;
    row.final_res_phi_fraction = result.final_res_phi_fraction;

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    row.source_region_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, body_length,
                                     body_height, body_width, core_shell, source_region);
    const Real exterior_solution_block_max =
        coreBlockMaxAbsNormOutsideRegion(particles, names.solution, positions, total_real_particles, body_length,
                                         body_height, body_width, core_shell, source_region);
    row.source_exterior_ratio = row.source_region_solution_block_max / (exterior_solution_block_max + TinyReal);

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names);
    row.global_div_a = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles);
    row.source_region_div_a = hostBodyDivAReductionMetrics(
        particles, names, positions, total_real_particles,
        [&](const Vecd &position) { return benchmark::insideBoxRegion(position, source_region); });
    return row;
}

inline void printInnerADivergencePenaltySweepRow(const char *test_name, const AphiInnerADivergencePenaltySweepRow &row)
{
    std::cout << test_name << " use_a_penalty=" << (row.use_a_divergence_penalty ? 1 : 0)
              << " a_penalty=" << row.a_divergence_penalty << " restart=" << row.restart_dimension
              << " max_outer=" << row.max_outer_iterations << " converged=" << (row.converged ? 1 : 0)
              << " outer=" << row.outer_iterations << " rel=" << row.final_relative_residual
              << " true_rel=" << row.final_true_relative_residual << " res_A_frac=" << row.final_res_a_fraction
              << " res_phi_frac=" << row.final_res_phi_fraction
              << " source_region_solution_block_max=" << row.source_region_solution_block_max
              << " source_exterior_ratio=" << row.source_exterior_ratio
              << " div_A_relative=" << row.global_div_a.div_a_relative
              << " source_div_A_relative=" << row.source_region_div_a.div_a_relative
              << " div_A_diag=pairwise"
              << " div_A_level=" << divAGaugeDiagnosticLevel(row.global_div_a.div_a_relative) << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_INNER_A_DIVERGENCE_PENALTY_SWEEP_HELPERS_H
