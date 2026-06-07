/**
 * Stage 10.8 Inner gate: solenoidal case at research eta_A vs plan §23 success criteria.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"

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
    const Real phi_gauge_penalty = 100.0;
    const Real tolerance = 1.0e-4;
    const Real max_observable_rel_change = 0.20;
    const Real solenoidal_current_amplitude = 5.0;
    const AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};
    const Real eta_baseline = 0.0;
    const Real eta_gate_values[] = {AphiADivergencePenaltyResearchDefaults::eta_a_min,
                                    AphiADivergencePenaltyResearchDefaults::eta_a_max};

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    const AphiInnerADivergencePenaltyGateRow baseline_row = runInnerSolenoidalGateRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_baseline, scale_metrics, tolerance, source_region, solenoidal_current_amplitude);
    printInnerADivergencePenaltyGateRow("test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic_baseline",
                                        baseline_row);

    size_t passed_gate_rows = 0;
    for (const Real eta_a : eta_gate_values)
    {
        const AphiInnerADivergencePenaltyGateRow row = runInnerSolenoidalGateRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, scale_metrics, tolerance, source_region, solenoidal_current_amplitude);
        printInnerADivergencePenaltyGateRow("test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic", row);

        const bool row_passed = innerSolenoidalGateRowPassed(
            row, tolerance, baseline_row.global_div_a.div_a_relative, baseline_row.global_joule_power,
            baseline_row.global_E_L2, max_observable_rel_change, eta_a > TinyReal);
        if (row_passed)
        {
            passed_gate_rows += 1;
        }
        std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic eta_A=" << eta_a
                  << " row_passed=" << (row_passed ? 1 : 0) << std::endl;
    }

    const bool baseline_ok = baseline_row.converged && baseline_row.final_true_relative_residual <= tolerance;
    const bool passed = baseline_ok && passed_gate_rows == (sizeof(eta_gate_values) / sizeof(eta_gate_values[0]));

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_gate_diagnostic baseline_ok=" << (baseline_ok ? 1 : 0)
              << " passed_gate_rows=" << passed_gate_rows << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
