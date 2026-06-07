/**
 * Stage 10.10-A/C: non-degenerate divergence-free MMS (Az2D, CrossSine3D) + exact E/J/Joule reference.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h"

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
    const Real max_block_linf_consistency = 0.15;
    const Real max_block_linf_invariance = 0.25;
    const Real min_rhs_norm = 1.0e-6;
    const Real max_joule_error_vs_exact = 0.25;
    const Real max_E_combined_error_vs_exact = 0.25;
    const Real max_B_error_vs_exact = 0.25;
    const Real max_J_combined_error_vs_exact = 0.25;
    const Real eta_values[] = {0.0, 0.1, 0.2, 0.3};
    const AphiDivFreeValidationFieldKind field_kinds[] = {AphiDivFreeValidationFieldKind::Az2D,
                                                        AphiDivFreeValidationFieldKind::CrossSine3D};

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    StdVec<AphiNonDegenerateDivFreeMMSRow> matrix_rows;
    matrix_rows.reserve(16);

    size_t consistency_passed = 0;
    size_t invariance_passed = 0;
    size_t expected_consistency = 0;
    size_t expected_invariance = 0;
    Real min_rhs_norm_seen = 1.0e30;
    const AphiNonDegenerateDivFreeMMSRow *az2d_eta0_row = nullptr;

    for (const AphiDivFreeValidationFieldKind field_kind : field_kinds)
    {
        for (const Real eta_a : eta_values)
        {
            expected_consistency += 1;
            const AphiNonDegenerateDivFreeMMSRow consistency_row = runNonDegenerateDivFreeMMSRow(
                ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty,
                eta_a, scale_metrics, tolerance, field_kind, true);
            matrix_rows.push_back(consistency_row);
            printNonDegenerateDivFreeMMSRow("test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms", consistency_row);
            min_rhs_norm_seen = std::min(min_rhs_norm_seen, consistency_row.rhs_norm);
            if (nonDegenerateDivFreeMMSConsistencyRowPassed(consistency_row, tolerance, max_block_linf_consistency,
                                                            min_rhs_norm))
            {
                consistency_passed += 1;
            }
            if (field_kind == AphiDivFreeValidationFieldKind::Az2D && eta_a <= TinyReal)
            {
                az2d_eta0_row = &matrix_rows.back();
            }

            if (eta_a <= TinyReal)
            {
                continue;
            }
            expected_invariance += 1;
            const AphiNonDegenerateDivFreeMMSRow invariance_row = runNonDegenerateDivFreeMMSRow(
                ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty,
                eta_a, scale_metrics, tolerance, field_kind, false);
            matrix_rows.push_back(invariance_row);
            printNonDegenerateDivFreeMMSRow("test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms", invariance_row);
            if (invariance_row.converged && invariance_row.rhs_norm >= min_rhs_norm &&
                std::isfinite(invariance_row.true_rel) && invariance_row.true_rel <= tolerance)
            {
                invariance_passed += 1;
            }
        }
    }

    printInnerEtaObservableMatrixSummary("test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms", matrix_rows);

    const bool exact_observable_ok =
        az2d_eta0_row != nullptr && az2d_eta0_row->joule_error_vs_exact <= max_joule_error_vs_exact &&
        az2d_eta0_row->E_combined_error_vs_exact <= max_E_combined_error_vs_exact &&
        az2d_eta0_row->b_error_vs_exact <= max_B_error_vs_exact &&
        az2d_eta0_row->J_combined_error_vs_exact <= max_J_combined_error_vs_exact;

    const bool passed = min_rhs_norm_seen >= min_rhs_norm && consistency_passed == expected_consistency &&
                        invariance_passed == expected_invariance && exact_observable_ok;

    std::cout << "test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms min_rhs_norm_seen=" << min_rhs_norm_seen
              << " consistency_passed=" << consistency_passed << "/" << expected_consistency
              << " invariance_passed=" << invariance_passed << "/" << expected_invariance
              << " az2d_eta0_joule_err_vs_exact=" << (az2d_eta0_row != nullptr ? az2d_eta0_row->joule_error_vs_exact : -1.0)
              << " az2d_eta0_E_combined_err_vs_exact="
              << (az2d_eta0_row != nullptr ? az2d_eta0_row->E_combined_error_vs_exact : -1.0)
              << " az2d_eta0_B_err_vs_exact=" << (az2d_eta0_row != nullptr ? az2d_eta0_row->b_error_vs_exact : -1.0)
              << " az2d_eta0_J_combined_err_vs_exact="
              << (az2d_eta0_row != nullptr ? az2d_eta0_row->J_combined_error_vs_exact : -1.0)
              << " exact_observable_ok=" << (exact_observable_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
