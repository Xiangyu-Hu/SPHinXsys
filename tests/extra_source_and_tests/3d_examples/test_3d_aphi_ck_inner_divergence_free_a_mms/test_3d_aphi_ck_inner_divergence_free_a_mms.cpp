/**
 * Stage 10.9-B: manufactured separable A-phi MMS + GMRES with optional A-penalty (consistency / invariance).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real tolerance = 1.0e-4;
    const Real max_block_linf_error = 0.15;
    const Real eta_values[] = {0.0, 0.1, 0.2, 0.3};

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    size_t consistency_passed = 0;
    size_t invariance_passed = 0;
    for (const Real eta_a : eta_values)
    {
        const AphiDivergenceFreeAMMSRow consistency_row = runManufacturedSeparableAphiMMSRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, scale_metrics, tolerance, true);
        printDivergenceFreeAMMSRow("test_3d_aphi_ck_inner_divergence_free_a_mms", consistency_row);
        if (divergenceFreeAMMSRowPassed(consistency_row, tolerance, max_block_linf_error))
        {
            consistency_passed += 1;
        }

        if (eta_a <= TinyReal)
        {
            continue;
        }
        const AphiDivergenceFreeAMMSRow invariance_row = runManufacturedSeparableAphiMMSRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, scale_metrics, tolerance, false);
        printDivergenceFreeAMMSRow("test_3d_aphi_ck_inner_divergence_free_a_mms", invariance_row);
        if (divergenceFreeAMMSInvarianceRowPassed(invariance_row, tolerance))
        {
            invariance_passed += 1;
        }
    }

    const size_t expected_consistency = sizeof(eta_values) / sizeof(eta_values[0]);
    const size_t expected_invariance = expected_consistency - 1;
    const bool passed = consistency_passed == expected_consistency && invariance_passed == expected_invariance;

    std::cout << "test_3d_aphi_ck_inner_divergence_free_a_mms manufactured_separable=1 consistency_passed="
              << consistency_passed << "/" << expected_consistency << " invariance_passed=" << invariance_passed << "/"
              << expected_invariance << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
