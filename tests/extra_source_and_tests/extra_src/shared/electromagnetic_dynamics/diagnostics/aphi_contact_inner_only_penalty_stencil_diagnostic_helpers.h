#ifndef APHI_CONTACT_INNER_ONLY_PENALTY_STENCIL_DIAGNOSTIC_HELPERS_H
#define APHI_CONTACT_INNER_ONLY_PENALTY_STENCIL_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_contact_interface_mms_consistency_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline const char *contactPenaltyStencilModeName(AphiContactADivergencePenaltyStencilMode mode)
{
    switch (mode)
    {
    case AphiContactADivergencePenaltyStencilMode::InnerContact:
        return "InnerContact";
    case AphiContactADivergencePenaltyStencilMode::InnerOnly:
        return "InnerOnly";
    default:
        return "Unknown";
    }
}

struct AphiContactPenaltyStencilComparisonRow
{
    AphiContactADivergencePenaltyStencilMode stencil = AphiContactADivergencePenaltyStencilMode::InnerContact;
    AphiContactInterfaceMmsConsistencyRow mms{};
    AphiContactDivFreeTwoBodyMmsRow gmres{};
};

struct AphiContactInnerOnlyPenaltyStencilDiagnosticSummary
{
    AphiContactInterfaceMmsConsistencyRow inner_only_eta0{};
    AphiContactPenaltyStencilComparisonRow inner_contact_eta01{};
    AphiContactPenaltyStencilComparisonRow inner_only_eta01{};
};

inline AphiContactPenaltyStencilComparisonRow runContactPenaltyStencilComparisonRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, Real x_interface, Real interface_band_half_width,
    AphiDivFreeValidationFieldKind field_kind, AphiContactADivergencePenaltyStencilMode penalty_stencil)
{
    AphiContactPenaltyStencilComparisonRow row;
    row.stencil = penalty_stencil;
    row.mms = runContactInterfaceMmsConsistencyRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_a, scale_metrics, x_interface, interface_band_half_width, field_kind, penalty_stencil);
    row.gmres = runContactDivFreeTwoBodyMmsRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_a, scale_metrics, tolerance, x_interface, interface_band_half_width, field_kind, 50, 100,
        penalty_stencil);
    return row;
}

inline void printContactPenaltyStencilComparisonRow(const char *test_name,
                                                    const AphiContactPenaltyStencilComparisonRow &row)
{
    std::cout << test_name << " stencil=" << contactPenaltyStencilModeName(row.stencil) << " eta_A=" << row.mms.eta_a
              << " mms_global_rel_l2=" << row.mms.global_relative_l2
              << " mms_core_safe_rel_l2=" << row.mms.stencil_safe_core.relative_l2
              << " mms_interface_rel_l2=" << row.mms.interface_band.relative_l2
              << " mms_boundary_rel_l2=" << row.mms.boundary.relative_l2
              << " gmres_converged=" << (row.gmres.converged ? 1 : 0)
              << " gmres_exact_consistency_defect=" << row.gmres.exact_consistency_defect
              << " gmres_global_true_rel=" << row.gmres.global_true_rel
              << " gmres_interface_band_rel=" << row.gmres.interface_band_rel << std::endl;
}

inline bool contactInnerOnlyPenaltyStencilEta0RegressionOk(const AphiContactInterfaceMmsConsistencyRow &row,
                                                             Real max_global_rel_l2, Real max_core_safe_rel_l2)
{
    return row.global_relative_l2 <= max_global_rel_l2 && row.stencil_safe_core.relative_l2 <= max_core_safe_rel_l2;
}

inline bool contactInnerOnlyPenaltyStencilResearchImproved(const AphiContactPenaltyStencilComparisonRow &baseline,
                                                          const AphiContactPenaltyStencilComparisonRow &inner_only,
                                                          Real min_inner_contact_worse_ratio = 1.0e3)
{
    if (inner_only.mms.eta_a <= TinyReal)
        return false;
    const Real global_ratio = baseline.mms.global_relative_l2 /
                              (inner_only.mms.global_relative_l2 + TinyReal);
    const Real gmres_defect_ratio = baseline.gmres.exact_consistency_defect /
                                    (inner_only.gmres.exact_consistency_defect + TinyReal);
    const bool global_mms_improved = global_ratio >= min_inner_contact_worse_ratio;
    const bool gmres_converged = inner_only.gmres.converged && std::isfinite(inner_only.gmres.global_true_rel);
    const bool gmres_defect_improved = gmres_defect_ratio >= min_inner_contact_worse_ratio;
    return global_mms_improved && gmres_converged && gmres_defect_improved;
}

inline bool contactInnerOnlyPenaltyStencilEta01MmsOk(const AphiContactPenaltyStencilComparisonRow &inner_only_eta01,
                                                     Real max_mms_global_rel_l2 = 1.0e-6,
                                                     Real max_gmres_global_true_rel = 1.0e-5)
{
    return inner_only_eta01.mms.eta_a > TinyReal && inner_only_eta01.gmres.converged &&
           std::isfinite(inner_only_eta01.gmres.global_true_rel) &&
           inner_only_eta01.mms.global_relative_l2 <= max_mms_global_rel_l2 &&
           inner_only_eta01.gmres.global_true_rel <= max_gmres_global_true_rel;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_INNER_ONLY_PENALTY_STENCIL_DIAGNOSTIC_HELPERS_H
