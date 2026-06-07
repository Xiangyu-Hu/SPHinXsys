#ifndef APHI_INNER_A_DIVERGENCE_PENALTY_GATE_HELPERS_H
#define APHI_INNER_A_DIVERGENCE_PENALTY_GATE_HELPERS_H

#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_eta_sweep_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiInnerADivergencePenaltyGateRow : public AphiInnerADivergencePenaltyEtaSweepRow
{
    Real global_joule_power = 0.0;
    Real global_E_L2 = 0.0;
    Real global_J_L2 = 0.0;
};

inline void execInnerJoulePostProcess(AphiLhsTestBody &test_body, const AphiVariableNames &names, Real omega,
                                      const AphiJouleHeatingFieldNames &joule_fields)
{
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);

    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    test_body.updateRelations();
}

inline AphiInnerADivergencePenaltyGateRow runInnerSolenoidalGateRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, const benchmark::AphiBoxRegion &source_region,
    Real solenoidal_current_amplitude, Real imag_to_real_ratio = 0.1, UnsignedInt restart_dimension = 50,
    UnsignedInt max_outer_iterations = 100)
{
    AphiInnerADivergencePenaltyGateRow row;
    row.eta_a = eta_a;
    row.restart_dimension = restart_dimension;
    row.max_outer_iterations = max_outer_iterations;
    row.use_a_divergence_penalty = eta_a > TinyReal;
    row.a_divergence_penalty = lambdaAFromEtaA(eta_a, scale_metrics);
    row.penalty_to_laplace_diag_ratio = eta_a;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    const AphiJouleHeatingFieldNames joule_fields;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignSolenoidalCurlCurrentRhsCK> assign_solenoidal_current(
        test_body.body, names.rhs, source_region, solenoidal_current_amplitude, imag_to_real_ratio);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = row.use_a_divergence_penalty;
    options.a_divergence_penalty = row.a_divergence_penalty;

    const AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    (void)register_joule_fields;
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

    execInnerJoulePostProcess(test_body, names, omega, joule_fields);

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

    const auto all_particles = [](const Vecd &) { return true; };
    const AphiRegionalElectromagneticObservables observables = hostRegionalElectromagneticObservablesFromJouleFields(
        particles, names, joule_fields, positions, total_real_particles, body_length, body_height, body_width,
        core_shell, all_particles, sigma);
    row.global_joule_power = observables.Joule_power;
    row.global_E_L2 = observables.E_L2;
    row.global_J_L2 = observables.J_L2;
    return row;
}

inline void printInnerADivergencePenaltyGateRow(const char *test_name, const AphiInnerADivergencePenaltyGateRow &row)
{
    printInnerADivergencePenaltyEtaSweepRow(test_name, row);
    std::cout << test_name << " gate_observables eta_A=" << row.eta_a << " joule_power=" << row.global_joule_power
              << " E_L2=" << row.global_E_L2 << " J_L2=" << row.global_J_L2 << std::endl;
}

inline bool innerSolenoidalGateRowPassed(const AphiInnerADivergencePenaltyGateRow &row, Real tolerance,
                                         Real baseline_div_a_relative, Real baseline_joule_power, Real baseline_E_L2,
                                         Real max_observable_rel_change, bool require_div_a_drop)
{
    if (!row.converged || row.final_true_relative_residual > tolerance)
    {
        return false;
    }
    if (require_div_a_drop && row.global_div_a.div_a_relative >= baseline_div_a_relative - TinyReal)
    {
        return false;
    }
    if (relativeMetricChange(baseline_joule_power, row.global_joule_power) > max_observable_rel_change)
    {
        return false;
    }
    if (relativeMetricChange(baseline_E_L2, row.global_E_L2) > max_observable_rel_change)
    {
        return false;
    }
    return true;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_INNER_A_DIVERGENCE_PENALTY_GATE_HELPERS_H
