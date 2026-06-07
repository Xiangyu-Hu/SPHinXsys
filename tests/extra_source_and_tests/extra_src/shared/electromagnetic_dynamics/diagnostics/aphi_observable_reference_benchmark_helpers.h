#ifndef APHI_OBSERVABLE_REFERENCE_BENCHMARK_HELPERS_H
#define APHI_OBSERVABLE_REFERENCE_BENCHMARK_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiObservableReferenceMetrics
{
    Real dp = 0.0;
    Real eta_a = 0.0;
    bool converged = false;
    Real true_rel = 0.0;
    Real div_a_relative = 0.0;
    Real joule_power = 0.0;
    Real E_L2 = 0.0;
    Real J_L2 = 0.0;
    Real block_linf_error = 0.0;
    bool is_mms = false;
};

struct AphiObservableReferenceComparison
{
    Real eta_a = 0.0;
    Real joule_error_vs_reference = 0.0;
    Real E_L2_error_vs_reference = 0.0;
    Real J_L2_error_vs_reference = 0.0;
    Real div_a_relative = 0.0;
    Real div_a_reduction_vs_reference = 0.0;
    bool converged = false;
    bool div_a_lower_than_reference = false;
};

inline AphiObservableReferenceMetrics solenoidalGateRowToObservableMetrics(const AphiInnerADivergencePenaltyGateRow &row,
                                                                           Real dp)
{
    AphiObservableReferenceMetrics metrics;
    metrics.dp = dp;
    metrics.eta_a = row.eta_a;
    metrics.converged = row.converged;
    metrics.true_rel = row.final_true_relative_residual;
    metrics.div_a_relative = row.global_div_a.div_a_relative;
    metrics.joule_power = row.global_joule_power;
    metrics.E_L2 = row.global_E_L2;
    metrics.J_L2 = row.global_J_L2;
    metrics.is_mms = false;
    return metrics;
}

inline AphiObservableReferenceMetrics runSolenoidalObservableReferenceRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, Real omega, Real phi_gauge_penalty, Real eta_a, const AphiCoreOperatorScaleMetrics &scale_metrics,
    Real tolerance, const benchmark::AphiBoxRegion &source_region, Real solenoidal_current_amplitude)
{
    const Real boundary_width = 3.0 * dp_0;
    const AphiInnerADivergencePenaltyGateRow row = runInnerSolenoidalGateRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_a, scale_metrics, tolerance, source_region, solenoidal_current_amplitude);
    return solenoidalGateRowToObservableMetrics(row, dp_0);
}

inline AphiObservableReferenceMetrics runManufacturedSeparableMMSObservableRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, Real omega, Real phi_gauge_penalty, Real eta_a, const AphiCoreOperatorScaleMetrics &scale_metrics,
    Real tolerance)
{
    const Real boundary_width = 3.0 * dp_0;
    const AphiDivergenceFreeAMMSRow mms_row = runManufacturedSeparableAphiMMSRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_a, scale_metrics, tolerance, true);

    AphiObservableReferenceMetrics metrics;
    metrics.dp = dp_0;
    metrics.eta_a = eta_a;
    metrics.converged = mms_row.converged;
    metrics.true_rel = mms_row.true_rel;
    metrics.div_a_relative = mms_row.div_a.div_a_relative;
    metrics.block_linf_error = mms_row.block_linf_error;
    metrics.joule_power = mms_row.joule_power;
    metrics.E_L2 = mms_row.E_L2;
    metrics.J_L2 = mms_row.J_L2;
    metrics.is_mms = true;
    return metrics;
}

inline AphiObservableReferenceComparison compareObservableReferenceMetrics(const AphiObservableReferenceMetrics &reference,
                                                                           const AphiObservableReferenceMetrics &candidate)
{
    AphiObservableReferenceComparison comparison;
    comparison.eta_a = candidate.eta_a;
    comparison.converged = candidate.converged;
    comparison.joule_error_vs_reference = relativeMetricChange(reference.joule_power, candidate.joule_power);
    comparison.E_L2_error_vs_reference = relativeMetricChange(reference.E_L2, candidate.E_L2);
    comparison.J_L2_error_vs_reference = relativeMetricChange(reference.J_L2, candidate.J_L2);
    comparison.div_a_relative = candidate.div_a_relative;
    comparison.div_a_reduction_vs_reference =
        reference.div_a_relative / (candidate.div_a_relative + TinyReal);
    comparison.div_a_lower_than_reference = candidate.div_a_relative < reference.div_a_relative - TinyReal;
    return comparison;
}

inline void printObservableReferenceMetrics(const char *test_name, const char *track, const AphiObservableReferenceMetrics &metrics)
{
    std::cout << test_name << " track=" << track << " kind=" << (metrics.is_mms ? "mms" : "solenoidal") << " dp="
              << metrics.dp << " eta_A=" << metrics.eta_a << " converged=" << (metrics.converged ? 1 : 0)
              << " true_rel=" << metrics.true_rel << " div_A_rel=" << metrics.div_a_relative
              << " joule_power=" << metrics.joule_power << " E_L2=" << metrics.E_L2 << " J_L2=" << metrics.J_L2;
    if (metrics.is_mms)
    {
        std::cout << " block_Linf_err=" << metrics.block_linf_error;
    }
    std::cout << std::endl;
}

inline void printObservableReferenceComparison(const char *test_name, const char *track,
                                               const AphiObservableReferenceComparison &comparison)
{
    std::cout << test_name << " track=" << track << " eta_A=" << comparison.eta_a << " converged="
              << (comparison.converged ? 1 : 0) << " joule_err_vs_ref=" << comparison.joule_error_vs_reference
              << " E_L2_err_vs_ref=" << comparison.E_L2_error_vs_reference
              << " J_L2_err_vs_ref=" << comparison.J_L2_error_vs_reference << " div_A_rel=" << comparison.div_a_relative
              << " div_A_reduction_vs_ref=" << comparison.div_a_reduction_vs_reference
              << " div_A_lower=" << (comparison.div_a_lower_than_reference ? 1 : 0) << std::endl;
}

inline bool observableReferenceCandidatePassed(const AphiObservableReferenceComparison &comparison, Real tolerance,
                                               Real max_observable_error_vs_reference, bool require_div_a_drop)
{
    if (!comparison.converged)
    {
        return false;
    }
    if (comparison.joule_error_vs_reference > max_observable_error_vs_reference)
    {
        return false;
    }
    if (comparison.E_L2_error_vs_reference > max_observable_error_vs_reference)
    {
        return false;
    }
    if (require_div_a_drop && !comparison.div_a_lower_than_reference)
    {
        return false;
    }
    (void)tolerance;
    return true;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_OBSERVABLE_REFERENCE_BENCHMARK_HELPERS_H
