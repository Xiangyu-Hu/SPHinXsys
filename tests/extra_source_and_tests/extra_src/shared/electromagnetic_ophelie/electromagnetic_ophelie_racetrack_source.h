#ifndef ELECTROMAGNETIC_OPHELIE_RACETRACK_SOURCE_H
#define ELECTROMAGNETIC_OPHELIE_RACETRACK_SOURCE_H

#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_team7_native_geometry.h"
#include "electromagnetic_ophelie_team7_probe.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

enum class OphelieCoilSourceModel
{
    VolumeETheta,
    FilamentRacetrack,
    VolumeRacetrack,
    SurfaceETheta,
    SurfaceRacetrack
};

inline const char *ophelieCoilSourceModelName(OphelieCoilSourceModel model)
{
    switch (model)
    {
    case OphelieCoilSourceModel::FilamentRacetrack:
        return "filament-racetrack";
    case OphelieCoilSourceModel::VolumeRacetrack:
        return "volume-racetrack";
    case OphelieCoilSourceModel::SurfaceETheta:
        return "surface-etheta";
    case OphelieCoilSourceModel::SurfaceRacetrack:
        return "surface-racetrack";
    default:
        return "volume-etheta";
    }
}

inline OphelieCoilSourceModel parseOphelieCoilSourceModel(const std::string &name)
{
    if (name == "filament-racetrack")
    {
        return OphelieCoilSourceModel::FilamentRacetrack;
    }
    if (name == "volume-racetrack")
    {
        return OphelieCoilSourceModel::VolumeRacetrack;
    }
    if (name == "surface-etheta")
    {
        return OphelieCoilSourceModel::SurfaceETheta;
    }
    if (name == "surface-racetrack")
    {
        return OphelieCoilSourceModel::SurfaceRacetrack;
    }
    return OphelieCoilSourceModel::VolumeETheta;
}

struct RacetrackSegment
{
    Vecd a = Vecd::Zero();
    Vecd b = Vecd::Zero();
    Vecd tangent = Vecd::Zero();
    Real length = 0.0;
};

struct FilamentQuadraturePoint
{
    Vecd position = Vecd::Zero();
    Vecd dl = Vecd::Zero();
};

struct OphelieRacetrackParams
{
    Real inset_mm = 0.0;
    Real z_mm = 99.0;
    Real ds_mm = 2.0;
};

struct Team7BzDiagnosticSummary
{
    OphelieCoilSourceModel source_model = OphelieCoilSourceModel::VolumeETheta;
    OphelieRacetrackParams racetrack;
    Team7BzCompareMetrics metrics;
    Real peak_ref_x_mm = 0.0;
    Real peak_sim_x_mm = 0.0;
    Real bz_at_126_sim_mT = 0.0;
    Real bz_at_126_ref_mT = 0.0;
    Real bz_at_198_sim_mT = 0.0;
    Real bz_at_198_ref_mT = 0.0;
    Real negative_lobe_rms_rel = 0.0;
    Real right_overshoot_mean_ratio = 0.0;
};

inline StdVec<RacetrackSegment> makeRectangularRacetrackFromStlBBox(const Team7NativeStlMeshBBoxMm &mesh_bbox,
                                                                    Real inset_mm, Real z_mm, Real mm_to_m = 1.0e-3)
{
    const Real x0 = mesh_bbox.coil_lower_[0] + inset_mm;
    const Real x1 = mesh_bbox.coil_upper_[0] - inset_mm;
    const Real y0 = mesh_bbox.coil_lower_[1] + inset_mm;
    const Real y1 = mesh_bbox.coil_upper_[1] - inset_mm;
    const Real z_m = z_mm * mm_to_m;
    const Vecd corners[4] = {Vecd(x0 * mm_to_m, y0 * mm_to_m, z_m), Vecd(x1 * mm_to_m, y0 * mm_to_m, z_m),
                             Vecd(x1 * mm_to_m, y1 * mm_to_m, z_m), Vecd(x0 * mm_to_m, y1 * mm_to_m, z_m)};
    StdVec<RacetrackSegment> segments;
    segments.reserve(4);
    for (int edge = 0; edge < 4; ++edge)
    {
        RacetrackSegment segment;
        segment.a = corners[edge];
        segment.b = corners[(edge + 1) % 4];
        const Vecd ab = segment.b - segment.a;
        segment.length = ab.norm();
        if (segment.length > TinyReal)
        {
            segment.tangent = ab / segment.length;
        }
        else
        {
            segment.tangent = Vecd::Zero();
        }
        segments.push_back(segment);
    }
    return segments;
}

inline void buildFilamentQuadrature(const StdVec<RacetrackSegment> &segments, Real ds_m,
                                    StdVec<FilamentQuadraturePoint> &quadrature)
{
    quadrature.clear();
    for (const RacetrackSegment &segment : segments)
    {
        if (segment.length <= TinyReal)
        {
            continue;
        }
        const int n_steps = std::max(1, static_cast<int>(std::ceil(segment.length / ds_m)));
        const Vecd dl = (segment.b - segment.a) / static_cast<Real>(n_steps);
        for (int step = 0; step < n_steps; ++step)
        {
            FilamentQuadraturePoint point;
            point.position = segment.a + (static_cast<Real>(step) + 0.5) * dl;
            point.dl = dl;
            quadrature.push_back(point);
        }
    }
}

inline Real racetrackPathLength(const StdVec<RacetrackSegment> &segments)
{
    Real length = 0.0;
    for (const RacetrackSegment &segment : segments)
    {
        length += segment.length;
    }
    return length;
}

inline Vecd nearestRacetrackTangent(const Vecd &position, const StdVec<RacetrackSegment> &segments)
{
    Real best_distance_sq = std::numeric_limits<Real>::infinity();
    Vecd best_tangent = Vecd::Zero();
    for (const RacetrackSegment &segment : segments)
    {
        const Vecd ab = segment.b - segment.a;
        const Real length_sq = ab.squaredNorm();
        if (length_sq <= TinyReal)
        {
            continue;
        }
        Real u = ab.dot(position - segment.a) / length_sq;
        u = std::max(Real(0.0), std::min(Real(1.0), u));
        const Vecd nearest = segment.a + u * ab;
        const Real distance_sq = (position - nearest).squaredNorm();
        if (distance_sq < best_distance_sq)
        {
            best_distance_sq = distance_sq;
            best_tangent = segment.tangent;
        }
    }
    return best_tangent;
}

inline void evaluateFilamentBiotSavartBzAtProbes(const StdVec<FilamentQuadraturePoint> &quadrature, Real total_current_a,
                                                 Real mu0, Real softening_length,
                                                 const StdVec<Team7BzProbePoint> &probes,
                                                 StdVec<Team7BzProbePoint> &results)
{
    results = probes;
    const Real coeff = mu0 / (4.0 * Pi);
    const Real eps2 = softening_length * softening_length;
    for (Team7BzProbePoint &probe : results)
    {
        Vecd b_sum = Vecd::Zero();
        for (const FilamentQuadraturePoint &source : quadrature)
        {
            const Vecd r = probe.position_m - source.position;
            const Real r2 = r.squaredNorm() + eps2;
            const Real inv_r3 = 1.0 / (std::sqrt(r2) * r2);
            const Vecd current_moment = total_current_a * source.dl;
            b_sum += coeff * current_moment.cross(r) * inv_r3;
        }
        probe.bz_sim_mT = b_sum[2] * 1000.0;
    }
}

inline Real probeBzAtXmm(const StdVec<Team7BzProbePoint> &probes, Real x_mm, bool simulated)
{
    for (const Team7BzProbePoint &probe : probes)
    {
        if (std::abs(probe.x_mm - x_mm) < 0.5)
        {
            return simulated ? probe.bz_sim_mT : probe.bz_ref_phase0_mT;
        }
    }
    return 0.0;
}

inline Team7BzDiagnosticSummary summarizeTeam7BzDiagnostic(const StdVec<Team7BzProbePoint> &probes,
                                                             const Team7BzCompareMetrics &metrics,
                                                             OphelieCoilSourceModel source_model,
                                                             const OphelieRacetrackParams &racetrack,
                                                             Real smoke_rms_threshold = 0.5)
{
    Team7BzDiagnosticSummary summary;
    summary.source_model = source_model;
    summary.racetrack = racetrack;
    summary.metrics = metrics;
    if (probes.empty())
    {
        return summary;
    }
    summary.peak_ref_x_mm = probes[metrics.peak_ref_index].x_mm;
    summary.peak_sim_x_mm = probes[metrics.peak_sim_index].x_mm;
    summary.bz_at_126_sim_mT = probeBzAtXmm(probes, 126.0, true);
    summary.bz_at_126_ref_mT = probeBzAtXmm(probes, 126.0, false);
    summary.bz_at_198_sim_mT = probeBzAtXmm(probes, 198.0, true);
    summary.bz_at_198_ref_mT = probeBzAtXmm(probes, 198.0, false);

    Real neg_ref_sq = 0.0;
    Real neg_err_sq = 0.0;
    Real right_ref_sum = 0.0;
    Real right_ratio_sum = 0.0;
    size_t right_count = 0;
    for (const Team7BzProbePoint &probe : probes)
    {
        if (probe.x_mm <= 72.0 && probe.bz_ref_phase0_mT < 0.0)
        {
            const Real err = probe.bz_sim_mT - probe.bz_ref_phase0_mT;
            neg_ref_sq += probe.bz_ref_phase0_mT * probe.bz_ref_phase0_mT;
            neg_err_sq += err * err;
        }
        if (probe.x_mm >= 162.0 && probe.x_mm <= 234.0)
        {
            right_ref_sum += std::abs(probe.bz_ref_phase0_mT);
            right_ratio_sum += (probe.bz_sim_mT - probe.bz_ref_phase0_mT) / (std::abs(probe.bz_ref_phase0_mT) + TinyReal);
            ++right_count;
        }
    }
    summary.negative_lobe_rms_rel = std::sqrt(neg_err_sq) / (std::sqrt(neg_ref_sq) + TinyReal);
    summary.right_overshoot_mean_ratio = right_count > 0 ? right_ratio_sum / static_cast<Real>(right_count) : 0.0;
    summary.metrics.passed = summary.metrics.n_probes > 0 && summary.metrics.rms_rel_error < smoke_rms_threshold;
    return summary;
}

inline void printTeam7BzDiagnosticSummary(const Team7BzDiagnosticSummary &summary)
{
    std::cout << "[ophelie] TEAM7 Bz diagnostic source_model=" << ophelieCoilSourceModelName(summary.source_model);
    if (summary.source_model == OphelieCoilSourceModel::FilamentRacetrack ||
        summary.source_model == OphelieCoilSourceModel::VolumeRacetrack)
    {
        std::cout << " inset_mm=" << summary.racetrack.inset_mm << " z_mm=" << summary.racetrack.z_mm
                  << " ds_mm=" << summary.racetrack.ds_mm;
    }
    std::cout << " RMS=" << summary.metrics.rms_rel_error << " peak_x_ref=" << summary.peak_ref_x_mm
              << " peak_x_sim=" << summary.peak_sim_x_mm << " Bz126_sim/ref=" << summary.bz_at_126_sim_mT << "/"
              << summary.bz_at_126_ref_mT << " Bz198_sim/ref=" << summary.bz_at_198_sim_mT << "/"
              << summary.bz_at_198_ref_mT << " neg_lobe_rms=" << summary.negative_lobe_rms_rel
              << " right_overshoot_mean=" << summary.right_overshoot_mean_ratio
              << " passed=" << (summary.metrics.passed ? 1 : 0) << std::endl;
}

inline Team7BzDiagnosticSummary evaluateFilamentRacetrackTeam7BzDiagnostic(
    const Team7NativeStlMeshBBoxMm &mesh_bbox, const OphelieRacetrackParams &racetrack, Real total_current_a,
    Real mu0, Real softening_length, const StdVec<Team7BzProbePoint> &reference_probes, Real smoke_rms_threshold)
{
    const StdVec<RacetrackSegment> segments =
        makeRectangularRacetrackFromStlBBox(mesh_bbox, racetrack.inset_mm, racetrack.z_mm);
    StdVec<FilamentQuadraturePoint> quadrature;
    buildFilamentQuadrature(segments, racetrack.ds_mm * 1.0e-3, quadrature);
    StdVec<Team7BzProbePoint> simulated_probes;
    evaluateFilamentBiotSavartBzAtProbes(quadrature, total_current_a, mu0, softening_length, reference_probes,
                                         simulated_probes);
    const Team7BzCompareMetrics metrics = compareTeam7BzPhase0(simulated_probes, smoke_rms_threshold);
    return summarizeTeam7BzDiagnostic(simulated_probes, metrics, OphelieCoilSourceModel::FilamentRacetrack, racetrack,
                                      smoke_rms_threshold);
}

inline void runFilamentRacetrackTeam7BzSweep(const Team7NativeStlMeshBBoxMm &mesh_bbox, Real total_current_a,
                                             Real mu0, Real softening_length,
                                             const StdVec<Team7BzProbePoint> &reference_probes, Real ds_mm,
                                             Real smoke_rms_threshold)
{
    const StdVec<Real> inset_values = {0.0, 10.0, 20.0, 30.0};
    const StdVec<Real> z_values = {49.0, 99.0, 149.0};
    std::cout << "[ophelie] TEAM7 filament-racetrack sweep ds_mm=" << ds_mm << std::endl;
    std::cout << "[ophelie]   source_model inset_mm z_mm ds_mm RMS peak_x_sim peak_x_ref Bz126_sim/ref Bz198_sim/ref "
                 "neg_lobe right_overshoot"
              << std::endl;
    for (const Real inset_mm : inset_values)
    {
        for (const Real z_mm : z_values)
        {
            OphelieRacetrackParams racetrack;
            racetrack.inset_mm = inset_mm;
            racetrack.z_mm = z_mm;
            racetrack.ds_mm = ds_mm;
            const Team7BzDiagnosticSummary summary = evaluateFilamentRacetrackTeam7BzDiagnostic(
                mesh_bbox, racetrack, total_current_a, mu0, softening_length, reference_probes, smoke_rms_threshold);
            std::cout << "[ophelie]   filament-racetrack " << inset_mm << " " << z_mm << " " << ds_mm << " "
                      << summary.metrics.rms_rel_error << " " << summary.peak_sim_x_mm << " " << summary.peak_ref_x_mm
                      << " " << summary.bz_at_126_sim_mT << "/" << summary.bz_at_126_ref_mT << " "
                      << summary.bz_at_198_sim_mT << "/" << summary.bz_at_198_ref_mT << " "
                      << summary.negative_lobe_rms_rel << " " << summary.right_overshoot_mean_ratio << std::endl;
        }
    }
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_RACETRACK_SOURCE_H
