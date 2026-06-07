#ifndef APHI_CONTACT_INNER_ONLY_GRADDIV_PC_CONSISTENCY_HELPERS_H
#define APHI_CONTACT_INNER_ONLY_GRADDIV_PC_CONSISTENCY_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_contact_graddiv_pc_diagnostic_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactInnerOnlyGradDivPcConsistencyResult
{
    AphiContactGradDivPcConsistencyMetrics inner_only_metrics{};
    AphiContactGradDivPcConsistencyMetrics inner_contact_metrics{};
    Real interface_pc_fd_gap_inner_contact_minus_inner_only = 0.0;
    bool contact_graddiv_skipped = false;
};

inline AphiContactInnerOnlyGradDivPcConsistencyResult runContactInnerOnlyGradDivPcConsistencyDiagnostic(
    int ac, char *av[], Real min_interface_pc_fd_gap = 1.0)
{
    AphiContactInnerOnlyGradDivPcConsistencyResult result;
    result.inner_only_metrics = runContactGradDivPcConsistencyMetrics(
        ac, av, AphiContactADivergencePenaltyStencilMode::InnerOnly);
    result.inner_contact_metrics = runContactGradDivPcConsistencyMetrics(
        ac, av, AphiContactADivergencePenaltyStencilMode::InnerContact);
    result.interface_pc_fd_gap_inner_contact_minus_inner_only =
        result.inner_contact_metrics.interface_pc_max_abs_diff - result.inner_only_metrics.interface_pc_max_abs_diff;
    result.contact_graddiv_skipped =
        result.inner_only_metrics.interface_pc_probed > 0 && result.inner_contact_metrics.interface_pc_probed > 0 &&
        result.interface_pc_fd_gap_inner_contact_minus_inner_only >= min_interface_pc_fd_gap;
    return result;
}

inline bool contactInnerOnlyGradDivPcConsistencyPassed(const AphiContactInnerOnlyGradDivPcConsistencyResult &result,
                                                       Real max_core_pc_diff = 5.0e-5)
{
    return result.inner_only_metrics.core_pc_probed > 0 &&
           result.inner_only_metrics.core_pc_max_abs_diff < max_core_pc_diff && result.contact_graddiv_skipped;
}

inline void printContactInnerOnlyGradDivPcConsistencyResult(const char *test_name,
                                                            const AphiContactInnerOnlyGradDivPcConsistencyResult &result)
{
    std::cout << test_name << " stencil=InnerOnly core_pc_probed=" << result.inner_only_metrics.core_pc_probed
              << " core_pc_max_abs_diff=" << result.inner_only_metrics.core_pc_max_abs_diff
              << " interface_pc_max_abs_diff=" << result.inner_only_metrics.interface_pc_max_abs_diff << std::endl;
    std::cout << test_name << " stencil=InnerContact core_pc_probed=" << result.inner_contact_metrics.core_pc_probed
              << " core_pc_max_abs_diff=" << result.inner_contact_metrics.core_pc_max_abs_diff
              << " interface_pc_max_abs_diff=" << result.inner_contact_metrics.interface_pc_max_abs_diff << std::endl;
    std::cout << test_name << " inner_only_pc_core_max_diff=" << result.inner_only_metrics.core_pc_max_abs_diff
              << " interface_pc_fd_gap_inner_contact_minus_inner_only="
              << result.interface_pc_fd_gap_inner_contact_minus_inner_only
              << " contact_graddiv_skipped=" << (result.contact_graddiv_skipped ? 1 : 0) << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_INNER_ONLY_GRADDIV_PC_CONSISTENCY_HELPERS_H
