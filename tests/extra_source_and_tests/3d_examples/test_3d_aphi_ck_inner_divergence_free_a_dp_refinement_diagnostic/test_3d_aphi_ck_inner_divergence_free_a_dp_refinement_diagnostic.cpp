/**
 * Stage 10.9-C: dp refinement for Linear2D pairwise divA + manufactured MMS consistency.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h"

#include <iostream>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real eta_a = 0.1;
    const Real tolerance = 1.0e-4;
    const Real max_linear2d_pairwise_div_a_rel = 0.01;
    const Real max_mms_block_linf_per_dp = 0.05;
    const Real dp_values[] = {0.15, 0.1, 0.075};

    std::vector<AphiDivergenceFreeADpRefinementRow> rows;
    rows.reserve(sizeof(dp_values) / sizeof(dp_values[0]));

    for (const Real dp_0 : dp_values)
    {
        const Real boundary_width = 3.0 * dp_0;
        const Real core_shell = 2.5 * dp_0;
        AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
        const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
            scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

        const AphiDivergenceFreeADpRefinementRow row = runDivergenceFreeADpRefinementRow(
            ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, eta_a,
            scale_metrics, tolerance);
        printDivergenceFreeADpRefinementRow("test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic", row);
        rows.push_back(row);
    }

    const bool passed = divergenceFreeADpRefinementPassed(rows, tolerance, max_linear2d_pairwise_div_a_rel,
                                                          max_mms_block_linf_per_dp);
    const bool mms_block_monotone = rows.back().mms_block_linf_error <= rows.front().mms_block_linf_error;
    std::cout << "test_3d_aphi_ck_inner_divergence_free_a_dp_refinement_diagnostic coarse_dp=" << rows.front().dp
              << " fine_dp=" << rows.back().dp << " coarse_mms_block_Linf=" << rows.front().mms_block_linf_error
              << " fine_mms_block_Linf=" << rows.back().mms_block_linf_error
              << " mms_block_monotone=" << (mms_block_monotone ? 1 : 0)
              << " fine_linear2d_pairwise_divA_rel=" << rows.back().linear2d_pairwise_div_a_rel
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
