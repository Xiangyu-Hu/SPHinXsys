/**
 * Stage 10.12: B=curl A dp refinement on exact divergence-free fields (B-corrected grad).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dp_refinement_diagnostic_helpers.h"

#include <iostream>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

std::vector<AphiCurlBDpRefinementRow> runFieldDpSweep(int ac, char *av[], const Real *dp_values, size_t dp_count,
                                                      Real body_length, Real body_height, Real body_width, Real sigma,
                                                      Real nu, Real omega, AphiDivFreeValidationFieldKind field_kind)
{
    std::vector<AphiCurlBDpRefinementRow> rows;
    rows.reserve(dp_count);
    for (size_t i = 0; i != dp_count; ++i)
    {
        const Real dp_0 = dp_values[i];
        const Real core_shell = 2.5 * dp_0;
        const AphiCurlBDpRefinementRow row =
            runCurlBDpRefinementRow(ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, omega,
                                    field_kind);
        printCurlBDpRefinementRow("test_3d_aphi_ck_curl_b_dp_refinement_diagnostic", row);
        rows.push_back(row);
    }
    return rows;
}

} // namespace

int main(int ac, char *av[])
{
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real max_linear2d_b_error = 0.01;
    const Real max_E_error = 0.05;
    const Real dp_values[] = {0.15, 0.10, 0.075, 0.05};

    const std::vector<AphiCurlBDpRefinementRow> linear2d_rows = runFieldDpSweep(
        ac, av, dp_values, sizeof(dp_values) / sizeof(dp_values[0]), body_length, body_height, body_width, sigma, nu,
        omega, AphiDivFreeValidationFieldKind::Linear2D);
    const std::vector<AphiCurlBDpRefinementRow> az2d_rows = runFieldDpSweep(
        ac, av, dp_values, sizeof(dp_values) / sizeof(dp_values[0]), body_length, body_height, body_width, sigma, nu,
        omega, AphiDivFreeValidationFieldKind::Az2D);
    const std::vector<AphiCurlBDpRefinementRow> crosssine_rows = runFieldDpSweep(
        ac, av, dp_values, sizeof(dp_values) / sizeof(dp_values[0]), body_length, body_height, body_width, sigma, nu,
        omega, AphiDivFreeValidationFieldKind::CrossSine3D);

    const Real max_fine_oscillatory_b_error = 0.05;
    const bool linear2d_ok = curlBDpRefinementLinear2DStable(linear2d_rows, max_linear2d_b_error);
    const bool az2d_converges = curlBDpRefinementOscillatoryConverges(az2d_rows, max_fine_oscillatory_b_error);
    const bool az2d_stage1013_ok = curlBDpRefinementAz2DStage1013Gate(az2d_rows);
    const bool crosssine_converges =
        curlBDpRefinementOscillatoryConverges(crosssine_rows, max_fine_oscillatory_b_error);

    const auto E_ok = [&](const std::vector<AphiCurlBDpRefinementRow> &rows) {
        for (const AphiCurlBDpRefinementRow &row : rows)
        {
            if (row.E_combined_error_vs_exact > max_E_error)
            {
                return false;
            }
        }
        return true;
    };

    const bool E_chain_ok = E_ok(linear2d_rows) && E_ok(az2d_rows) && E_ok(crosssine_rows);
    const bool passed = linear2d_ok && az2d_converges && az2d_stage1013_ok && crosssine_converges && E_chain_ok;

    const AphiCurlBDpRefinementRow *const az2d_dp10 = findCurlBDpRefinementRowByDp(az2d_rows, 0.10);
    const AphiCurlBDpRefinementRow *const az2d_dp075 = findCurlBDpRefinementRowByDp(az2d_rows, 0.075);
    const AphiCurlBDpRefinementRow *const az2d_dp05 = findCurlBDpRefinementRowByDp(az2d_rows, 0.05);

    std::cout << "test_3d_aphi_ck_curl_b_dp_refinement_diagnostic linear2d_ok=" << (linear2d_ok ? 1 : 0)
              << " az2d_B_converges=" << (az2d_converges ? 1 : 0)
              << " az2d_stage1013_ok=" << (az2d_stage1013_ok ? 1 : 0)
              << " crosssine_B_converges=" << (crosssine_converges ? 1 : 0) << " E_chain_ok=" << (E_chain_ok ? 1 : 0)
              << " az2d_B_dp0.10=" << (az2d_dp10 != nullptr ? az2d_dp10->b_error_vs_exact : -1.0)
              << " az2d_B_dp0.075=" << (az2d_dp075 != nullptr ? az2d_dp075->b_error_vs_exact : -1.0)
              << " az2d_B_dp0.05=" << (az2d_dp05 != nullptr ? az2d_dp05->b_error_vs_exact : -1.0)
              << " crosssine_coarse_B=" << crosssine_rows.front().b_error_vs_exact
              << " crosssine_medium_B=" << crosssine_rows[crosssine_rows.size() >= 3 ? 2 : 1].b_error_vs_exact
              << " crosssine_fine_B=" << crosssine_rows.back().b_error_vs_exact << " passed=" << (passed ? 1 : 0)
              << std::endl;
    (void)nu;
    return passed ? 0 : 1;
}
