/**
 * Stage 10.9-D / Case C: E/J/Joule observable reference benchmark.
 * Solenoidal: fine-dp eta=0 self-reference; MMS: same-dp eta=0 consistency reference.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_observable_reference_benchmark_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
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
    const Real tolerance = 1.0e-4;
    const Real solenoidal_current_amplitude = 5.0;
    const AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};

    const Real reference_dp = 0.075;
    const Real candidate_dp = 0.1;
    const Real reference_core_shell = 2.0 * reference_dp;
    const Real candidate_core_shell = 2.0 * candidate_dp;
    const Real mms_core_shell = 2.5 * candidate_dp;

    const Real max_solenoidal_observable_error_vs_cand_eta0 = 0.25;
    const Real max_mms_observable_error_vs_eta0 = 0.20;
    const Real eta_candidates[] = {AphiADivergencePenaltyResearchDefaults::eta_a_min, 0.2,
                                   AphiADivergencePenaltyResearchDefaults::eta_a_max};

    AphiLhsTestBody sol_ref_scale_body(reference_dp, body_length, body_height, body_width, 3.0 * reference_dp, ac, av);
    const AphiCoreOperatorScaleMetrics sol_ref_scale = hostCoreOperatorScaleMetrics(
        sol_ref_scale_body, body_length, body_height, body_width, reference_core_shell, sigma, nu, omega,
        phi_gauge_penalty, ac, av);
    const AphiObservableReferenceMetrics solenoidal_fine_reference = runSolenoidalObservableReferenceRow(
        ac, av, reference_dp, body_length, body_height, body_width, reference_core_shell, sigma, nu, omega,
        phi_gauge_penalty, 0.0, sol_ref_scale, tolerance, source_region, solenoidal_current_amplitude);
    printObservableReferenceMetrics("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                  "solenoidal_fine_ref", solenoidal_fine_reference);

    AphiLhsTestBody sol_cand_scale_body(candidate_dp, body_length, body_height, body_width, 3.0 * candidate_dp, ac, av);
    const AphiCoreOperatorScaleMetrics sol_cand_scale = hostCoreOperatorScaleMetrics(
        sol_cand_scale_body, body_length, body_height, body_width, candidate_core_shell, sigma, nu, omega,
        phi_gauge_penalty, ac, av);
    const AphiObservableReferenceMetrics solenoidal_cand_baseline = runSolenoidalObservableReferenceRow(
        ac, av, candidate_dp, body_length, body_height, body_width, candidate_core_shell, sigma, nu, omega,
        phi_gauge_penalty, 0.0, sol_cand_scale, tolerance, source_region, solenoidal_current_amplitude);
    printObservableReferenceMetrics("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                  "solenoidal_cand_baseline", solenoidal_cand_baseline);

    const AphiObservableReferenceComparison sol_mesh_gap =
        compareObservableReferenceMetrics(solenoidal_fine_reference, solenoidal_cand_baseline);
    printObservableReferenceComparison("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                       "solenoidal_mesh_gap_fine_vs_cand_eta0", sol_mesh_gap);

    size_t solenoidal_passed = 0;
    for (const Real eta_a : eta_candidates)
    {
        const AphiObservableReferenceMetrics candidate = runSolenoidalObservableReferenceRow(
            ac, av, candidate_dp, body_length, body_height, body_width, candidate_core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, sol_cand_scale, tolerance, source_region, solenoidal_current_amplitude);
        printObservableReferenceMetrics("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                        "solenoidal_candidate", candidate);

        const AphiObservableReferenceComparison vs_fine_ref =
            compareObservableReferenceMetrics(solenoidal_fine_reference, candidate);
        printObservableReferenceComparison("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                           "solenoidal_vs_fine_ref", vs_fine_ref);

        const AphiObservableReferenceComparison vs_cand_baseline =
            compareObservableReferenceMetrics(solenoidal_cand_baseline, candidate);
        printObservableReferenceComparison("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                           "solenoidal_vs_cand_eta0", vs_cand_baseline);

        if (observableReferenceCandidatePassed(vs_cand_baseline, tolerance,
                                               max_solenoidal_observable_error_vs_cand_eta0, true))
        {
            solenoidal_passed += 1;
        }
    }

    AphiLhsTestBody mms_scale_body(candidate_dp, body_length, body_height, body_width, 3.0 * candidate_dp, ac, av);
    const AphiCoreOperatorScaleMetrics mms_scale = hostCoreOperatorScaleMetrics(
        mms_scale_body, body_length, body_height, body_width, mms_core_shell, sigma, nu, omega, phi_gauge_penalty, ac,
        av);
    const AphiObservableReferenceMetrics mms_reference = runManufacturedSeparableMMSObservableRow(
        ac, av, candidate_dp, body_length, body_height, body_width, mms_core_shell, sigma, nu, omega,
        phi_gauge_penalty, 0.0, mms_scale, tolerance);
    printObservableReferenceMetrics("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                    "mms_eta0_ref", mms_reference);

    size_t mms_passed = 0;
    for (const Real eta_a : eta_candidates)
    {
        const AphiObservableReferenceMetrics candidate = runManufacturedSeparableMMSObservableRow(
            ac, av, candidate_dp, body_length, body_height, body_width, mms_core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, mms_scale, tolerance);
        printObservableReferenceMetrics("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                        "mms_candidate", candidate);

        const AphiObservableReferenceComparison vs_mms_ref = compareObservableReferenceMetrics(mms_reference, candidate);
        printObservableReferenceComparison("test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic",
                                           "mms_vs_eta0_ref", vs_mms_ref);

        if (observableReferenceCandidatePassed(vs_mms_ref, tolerance, max_mms_observable_error_vs_eta0, false))
        {
            mms_passed += 1;
        }
    }

    const size_t expected_rows = sizeof(eta_candidates) / sizeof(eta_candidates[0]);
    const bool reference_ok = solenoidal_fine_reference.converged && mms_reference.converged &&
                              solenoidal_cand_baseline.converged;
    const bool passed = reference_ok && solenoidal_passed == expected_rows && mms_passed == expected_rows;

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_observable_reference_diagnostic reference_ok="
              << (reference_ok ? 1 : 0) << " solenoidal_passed=" << solenoidal_passed << "/" << expected_rows
              << " mms_passed=" << mms_passed << "/" << expected_rows << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
