#ifndef APHI_TEAM7_CANONICAL_CASE_CK_H
#define APHI_TEAM7_CANONICAL_CASE_CK_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace benchmark
{

/**
 * Canonical TEAM7-like case specification (simplified fraction layout on physical box).
 *
 * Reference strategy (Stage 10B):
 * - Primary: fine-mesh self-reference at reference_dp (no external FEM in repo yet).
 * - Recorded observables from reference_dp run serve as regression anchors.
 * - External FEM/MFEM comparison: deferred until reference data is imported.
 *
 * Geometry note: regions use unit-box fractions mapped to TEAM7 physical dimensions
 * (1.2 x 1.0 x 0.3 m). This is a TEAM7-inspired benchmark, not a full CAD replica.
 */
struct AphiTeam7CanonicalCaseSpec
{
    static constexpr Real body_length = AphiTeam7PhysicalDimensions::length;
    static constexpr Real body_height = AphiTeam7PhysicalDimensions::height;
    static constexpr Real body_width = AphiTeam7PhysicalDimensions::width;

    static constexpr Real omega = 1.25;
    static constexpr Real phi_gauge_penalty = 100.0;
    static constexpr Real impressed_current_amplitude = 8.0;

    static constexpr Real tolerance = 5.0e-4;
    static constexpr UnsignedInt restart_dimension = 50;
    static constexpr UnsignedInt max_outer_iterations = 150;

    /** Fine-mesh self-reference resolution for quantitative comparison. */
    static constexpr Real reference_dp = 0.075;
    /** Engineering candidate resolution (matches 9E/10A baseline). */
    static constexpr Real candidate_dp = 0.1;

    static constexpr UnsignedInt centerline_bins = 16;
    static constexpr Real centerline_yz_band_factor = 1.5;

    /** Regression anchors from reference_dp run (2026-05-21, tol=5e-4, m=50 coupled PC). */
    static constexpr Real recorded_conductor_joule_power = 0.0931;
    static constexpr Real recorded_conductor_solution_block_max = 3.07;
    static constexpr Real recorded_coil_joule_power = 0.0;
    static constexpr Real recorded_total_joule_power = 0.136;

    /** Max relative observable change allowed for candidate_dp vs reference_dp (Joule + field). */
    static constexpr Real max_candidate_vs_reference_rel = 0.20;
    /** Centerline profile: informational only until mesh is finer / bins better populated. */
    static constexpr Real max_centerline_profile_rel = 2.0;
    /** Max relative drift allowed vs recorded regression anchors @ reference_dp. */
    static constexpr Real max_reference_vs_recorded_rel = 0.05;
};

} // namespace benchmark
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_CANONICAL_CASE_CK_H
