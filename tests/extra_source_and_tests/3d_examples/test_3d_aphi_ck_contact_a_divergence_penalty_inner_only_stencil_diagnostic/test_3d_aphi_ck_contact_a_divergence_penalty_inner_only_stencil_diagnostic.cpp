/**
 * Stage 10.12-C4: compare InnerContact vs InnerOnly penalty divA/grad(divA) stencil on two-body Contact.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_inner_only_penalty_stencil_diagnostic_helpers.h"

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
    const Real eta_a_research = 0.1;
    const Real max_global_rel_l2_eta0 = 1.0e-4;
    const Real max_core_safe_rel_l2_eta0 = 1.0e-4;
    const Real max_inner_only_eta01_mms_global = 1.0e-6;
    const Real max_inner_only_eta01_gmres_true_rel = 1.0e-5;
    const Real min_inner_contact_worse_ratio = 1.0e3;
    const AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    const char *test_name = "test_3d_aphi_ck_contact_a_divergence_penalty_inner_only_stencil_diagnostic";

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    const AphiContactInterfaceMmsConsistencyRow inner_only_eta0 = runContactInterfaceMmsConsistencyRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, 0.0, scale_metrics, x_interface, interface_band_half_width, field_kind,
        AphiContactADivergencePenaltyStencilMode::InnerOnly);
    std::cout << test_name << " stencil=InnerOnly eta_A=0 global_rel_l2=" << inner_only_eta0.global_relative_l2
              << " core_safe_rel_l2=" << inner_only_eta0.stencil_safe_core.relative_l2 << std::endl;

    const AphiContactPenaltyStencilComparisonRow inner_contact_eta01 = runContactPenaltyStencilComparisonRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_a_research, scale_metrics, tolerance, x_interface, interface_band_half_width, field_kind,
        AphiContactADivergencePenaltyStencilMode::InnerContact);
    const AphiContactPenaltyStencilComparisonRow inner_only_eta01 = runContactPenaltyStencilComparisonRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_a_research, scale_metrics, tolerance, x_interface, interface_band_half_width, field_kind,
        AphiContactADivergencePenaltyStencilMode::InnerOnly);

    printContactPenaltyStencilComparisonRow(test_name, inner_contact_eta01);
    printContactPenaltyStencilComparisonRow(test_name, inner_only_eta01);

    const bool eta0_regression_ok =
        contactInnerOnlyPenaltyStencilEta0RegressionOk(inner_only_eta0, max_global_rel_l2_eta0, max_core_safe_rel_l2_eta0);
    const bool comparison_reported = std::isfinite(inner_contact_eta01.mms.global_relative_l2) &&
                                     std::isfinite(inner_only_eta01.mms.global_relative_l2) &&
                                     std::isfinite(inner_contact_eta01.gmres.exact_consistency_defect) &&
                                     std::isfinite(inner_only_eta01.gmres.exact_consistency_defect);
    const bool research_improved = contactInnerOnlyPenaltyStencilResearchImproved(
        inner_contact_eta01, inner_only_eta01, min_inner_contact_worse_ratio);
    const bool inner_only_eta01_ok =
        contactInnerOnlyPenaltyStencilEta01MmsOk(inner_only_eta01, max_inner_only_eta01_mms_global,
                                                 max_inner_only_eta01_gmres_true_rel);
    const Real inner_contact_vs_inner_only_global_ratio =
        inner_contact_eta01.mms.global_relative_l2 / (inner_only_eta01.mms.global_relative_l2 + TinyReal);
    const bool passed =
        eta0_regression_ok && comparison_reported && research_improved && inner_only_eta01_ok;

    std::cout << test_name << " eta0_regression_ok=" << (eta0_regression_ok ? 1 : 0)
              << " comparison_reported=" << (comparison_reported ? 1 : 0)
              << " research_improved=" << (research_improved ? 1 : 0)
              << " inner_only_eta01_ok=" << (inner_only_eta01_ok ? 1 : 0)
              << " inner_contact_vs_inner_only_global_ratio=" << inner_contact_vs_inner_only_global_ratio
              << " passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
