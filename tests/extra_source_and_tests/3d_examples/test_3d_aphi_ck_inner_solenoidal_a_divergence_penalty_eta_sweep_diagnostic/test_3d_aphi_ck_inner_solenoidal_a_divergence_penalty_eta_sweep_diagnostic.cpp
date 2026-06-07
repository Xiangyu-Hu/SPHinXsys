/**
 * Stage 10.8: relative eta_A sweep for Inner solenoidal source (pairwise divA diagnostic).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_eta_sweep_helpers.h"

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
    const Real solenoidal_current_amplitude = 5.0;
    const AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};
    const Real eta_values[] = {0.0, 0.01, 0.03, 0.1, 0.3, 1.0};

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    std::cout << "test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_eta_sweep_diagnostic begin"
              << " median_laplace_a_diag=" << scale_metrics.median_abs_laplace_a_diag
              << " median_graddiv_block_norm=" << scale_metrics.median_graddiv_block_norm
              << " core_particles=" << scale_metrics.core_particles << std::endl;

    size_t converged_rows = 0;
    for (const Real eta_a : eta_values)
    {
        const AphiInnerADivergencePenaltyEtaSweepRow row = runInnerSolenoidalEtaSweepRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, scale_metrics, tolerance, source_region, solenoidal_current_amplitude);
        printInnerADivergencePenaltyEtaSweepRow(
            "test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_eta_sweep_diagnostic", row);
        if (row.converged)
        {
            converged_rows += 1;
        }
    }

    std::cout << "test_3d_aphi_ck_inner_solenoidal_a_divergence_penalty_eta_sweep_diagnostic completed_rows="
              << (sizeof(eta_values) / sizeof(eta_values[0])) << " converged_rows=" << converged_rows
              << " passed=1 (diagnostic: converged_rows not gating)" << std::endl;
    return 0;
}
