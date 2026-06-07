/**
 * Stage 10.12-C3+: two-body Contact A-penalty MMS (InnerOnly penalty stencil).
 * Strict gate: eta_A in {0, 0.1} (research primary). eta_A=0.2 informational.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_a_divergence_penalty_two_body_mms_helpers.h"

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
    const Real core_shell = dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real x_interface = 0.5;
    const Real interface_band_half_width = 2.0 * dp_0;
    const Real tolerance = 1.0e-4;
    const Real max_block_linf_error = 0.15;
    const Real max_exact_consistency_defect = 1.0e-4;
    const Real max_joule_error_vs_exact = 0.25;
    const Real max_E_combined_error_vs_exact = 0.25;
    const Real max_J_combined_error_vs_exact = 0.25;
    const Real max_upper_eta_exact_consistency_defect = 0.05;
    const Real eta_values[] = {0.0, 0.1, 0.2};
    const AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    const AphiContactADivergencePenaltyStencilMode penalty_stencil =
        AphiContactADivergencePenaltyStencilMode::InnerOnly;

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    size_t primary_rows_passed = 0;
    size_t expected_primary_rows = 0;
    size_t upper_eta_reported = 0;

    for (const Real eta_a : eta_values)
    {
        const AphiContactDivFreeTwoBodyMmsRow row = runContactDivFreeTwoBodyMmsRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, scale_metrics, tolerance, x_interface, interface_band_half_width, field_kind, 50,
            100, penalty_stencil);
        printContactDivFreeTwoBodyMmsRow("test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms", row);

        if (eta_a <= AphiADivergencePenaltyResearchDefaults::primary_eta_a + TinyReal)
        {
            expected_primary_rows += 1;
            if (contactDivFreeTwoBodyMmsRowStrictPassed(row, tolerance, max_block_linf_error,
                                                        max_exact_consistency_defect, max_joule_error_vs_exact,
                                                        max_E_combined_error_vs_exact, max_J_combined_error_vs_exact))
            {
                primary_rows_passed += 1;
            }
            continue;
        }

        if (std::isfinite(row.exact_consistency_defect) &&
            row.exact_consistency_defect <= max_upper_eta_exact_consistency_defect)
        {
            upper_eta_reported += 1;
        }
    }

    const bool primary_passed = primary_rows_passed == expected_primary_rows;
    const bool upper_eta_reported_ok = upper_eta_reported == 1;
    const bool passed = primary_passed && upper_eta_reported_ok;

    std::cout << "test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms penalty_stencil=InnerOnly"
              << " primary_rows_passed=" << primary_rows_passed << "/" << expected_primary_rows
              << " upper_eta_reported=" << upper_eta_reported << "/1 passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
