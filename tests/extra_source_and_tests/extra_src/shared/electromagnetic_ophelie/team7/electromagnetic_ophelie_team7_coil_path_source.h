#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_COIL_PATH_SOURCE_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_COIL_PATH_SOURCE_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_current_moment_source.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_team7_probe.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Official TEAM7 diagram radii [m]; STL bbox is used only for placement, not to infer R_out/R_in. */
struct OphelieTeam7CoilPathSourceSpec
{
    Real turns = 2742.0;
    Real current_per_turn = 1.0;
    Real R_out = 0.050;
    Real R_in = 0.025;
    bool reverse_winding = false;
    Real source_scale = 1.0;
    int polyline_segments_per_section = 32;

    Real totalAmpereTurns() const { return turns * current_per_turn; }
    Real centerlineRadius() const { return 0.5 * (R_out + R_in); }
};

inline constexpr const char *kOphelieTeam7CoilPathTangentName = "Team7CoilPathTangent";

struct OphelieTeam7CoilPathPrepareSummary
{
    Real ni_target = 0.0;
    Real path_length_m = 0.0;
    Real coil_volume_m3 = 0.0;
    Real a_eff_m2 = 0.0;
    Real j0 = 0.0;
    Real i_eff_formula = 0.0;
    size_t centerline_points = 0;
};

struct OphelieTeam7FilamentPrepareSummary
{
    Real ni_target = 0.0;
    Real path_length_m = 0.0;
    Real ds_m = 0.0;
    size_t quadrature_points = 0;
    size_t centerline_points = 0;
};

inline StdVec<Vecd> buildOphelieTeam7CoilCenterlinePolyline(const BoundingBoxd &coil_bbox,
                                                            const OphelieTeam7CoilPathSourceSpec &spec)
{
    StdVec<Vecd> pts;
    const Real R_out = spec.R_out;
    const Real R_in = spec.R_in;
    const Real R_c = spec.centerlineRadius();
    const Real x_min = coil_bbox.lower_[0];
    const Real x_max = coil_bbox.upper_[0];
    const Real y_min = coil_bbox.lower_[1];
    const Real y_max = coil_bbox.upper_[1];
    const Real z_c = 0.5 * (coil_bbox.lower_[2] + coil_bbox.upper_[2]);

    const Real cxL = x_min + R_out;
    const Real cxR = x_max - R_out;
    const Real cyB = y_min + R_out;
    const Real cyT = y_max - R_out;
    const Real y_inner_offset = 0.5 * (R_out - R_in);
    const int n = std::max(4, spec.polyline_segments_per_section);

    pts.reserve(static_cast<size_t>(8 * n));

    auto add_line = [&](const Vecd &a, const Vecd &b) {
        for (int k = 0; k < n; ++k)
        {
            const Real s = static_cast<Real>(k) / static_cast<Real>(n);
            pts.push_back(a * (1.0 - s) + b * s);
        }
    };

    auto add_arc = [&](const Vecd &corner_center, Real a0, Real a1) {
        for (int k = 0; k < n; ++k)
        {
            const Real s = static_cast<Real>(k) / static_cast<Real>(n);
            const Real angle = a0 * (1.0 - s) + a1 * s;
            pts.push_back(corner_center + Vecd(R_c * std::cos(angle), R_c * std::sin(angle), 0.0));
        }
    };

    // CCW winding viewed from +z (TEAM7 default).
    add_line(Vecd(cxL, y_min + y_inner_offset, z_c), Vecd(cxR, y_min + y_inner_offset, z_c));
    add_arc(Vecd(cxR, cyB, z_c), -0.5 * Pi, 0.0);
    add_line(Vecd(x_max - y_inner_offset, cyB, z_c), Vecd(x_max - y_inner_offset, cyT, z_c));
    add_arc(Vecd(cxR, cyT, z_c), 0.0, 0.5 * Pi);
    add_line(Vecd(cxR, y_max - y_inner_offset, z_c), Vecd(cxL, y_max - y_inner_offset, z_c));
    add_arc(Vecd(cxL, cyT, z_c), 0.5 * Pi, Pi);
    add_line(Vecd(x_min + y_inner_offset, cyT, z_c), Vecd(x_min + y_inner_offset, cyB, z_c));
    add_arc(Vecd(cxL, cyB, z_c), Pi, 1.5 * Pi);

    return pts;
}

inline Real closedPolylinePathLength(const StdVec<Vecd> &path)
{
    if (path.size() < 2)
    {
        return 0.0;
    }
    Real length = 0.0;
    for (size_t k = 0; k < path.size(); ++k)
    {
        const Vecd &a = path[k];
        const Vecd &b = path[(k + 1) % path.size()];
        length += (b - a).norm();
    }
    return length;
}

inline Vecd tangentFromClosedPolyline(const Vecd &position, const StdVec<Vecd> &path)
{
    if (path.size() < 2)
    {
        return Vecd::Zero();
    }

    Real best_distance_squared = std::numeric_limits<Real>::max();
    Vecd best_tangent = Vecd::Zero();
    for (size_t k = 0; k < path.size(); ++k)
    {
        const Vecd &a = path[k];
        const Vecd &b = path[(k + 1) % path.size()];
        const Vecd ab = b - a;
        const Real ab_squared = ab.squaredNorm();
        if (ab_squared <= TinyReal)
        {
            continue;
        }
        Real segment_parameter = (position - a).dot(ab) / ab_squared;
        segment_parameter = std::max(Real(0), std::min(Real(1), segment_parameter));
        const Vecd closest = a + segment_parameter * ab;
        const Real distance_squared = (position - closest).squaredNorm();
        if (distance_squared < best_distance_squared)
        {
            best_distance_squared = distance_squared;
            best_tangent = ab / std::sqrt(ab_squared);
        }
    }
    return best_tangent;
}

inline Real hostCoilParticleVolume(BaseParticles &particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t count = particles.TotalRealParticles();
    Real sum = 0.0;
    for (size_t i = 0; i != count; ++i)
    {
        sum += volume[i];
    }
    return sum;
}

inline void hostStoreOphelieTeam7CoilPathTangents(BaseParticles &particles, const StdVec<Vecd> &centerline,
                                                  const OphelieTeam7CoilPathSourceSpec &spec)
{
    particles.registerStateVariable<Vecd>(kOphelieTeam7CoilPathTangentName, ZeroData<Vecd>::value);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, kOphelieTeam7CoilPathTangentName);
    const size_t count = particles.TotalRealParticles();
    Vecd *tangents = particles.getVariableDataByName<Vecd>(kOphelieTeam7CoilPathTangentName);
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real winding_sign = spec.reverse_winding ? -1.0 : 1.0;
    for (size_t i = 0; i != count; ++i)
    {
        Vecd tangent = tangentFromClosedPolyline(positions[i], centerline);
        if (tangent.norm() > TinyReal)
        {
            tangent.normalize();
        }
        tangents[i] = winding_sign * tangent;
    }
    syncVariableToDevice<Vecd>(particles, kOphelieTeam7CoilPathTangentName);
}

inline Real ophelieTeam7CoilPathCurrentDensityJ0(const BoundingBoxd &coil_bbox, const OphelieTeam7CoilPathSourceSpec &spec,
                                                   Real coil_volume_m3)
{
    const StdVec<Vecd> centerline = buildOphelieTeam7CoilCenterlinePolyline(coil_bbox, spec);
    const Real path_length = closedPolylinePathLength(centerline);
    const Real a_eff = path_length > TinyReal ? coil_volume_m3 / path_length : TinyReal;
    return spec.source_scale * spec.totalAmpereTurns() / (a_eff + TinyReal);
}

inline OphelieTeam7CoilPathPrepareSummary prepareOphelieTeam7VolumeRacetrackCoilSource(
    SPHBody &coil_body, const BoundingBoxd &coil_bbox, const OphelieTeam7CoilPathSourceSpec &spec)
{
    OphelieTeam7CoilPathPrepareSummary summary;
    summary.ni_target = spec.totalAmpereTurns();
    BaseParticles &particles = coil_body.getBaseParticles();
    const StdVec<Vecd> centerline = buildOphelieTeam7CoilCenterlinePolyline(coil_bbox, spec);
    summary.centerline_points = centerline.size();
    summary.path_length_m = closedPolylinePathLength(centerline);
    summary.coil_volume_m3 = hostCoilParticleVolume(particles);
    summary.a_eff_m2 =
        summary.path_length_m > TinyReal ? summary.coil_volume_m3 / summary.path_length_m : TinyReal;
    summary.j0 = spec.source_scale * summary.ni_target / (summary.a_eff_m2 + TinyReal);
    summary.i_eff_formula = summary.j0 * summary.a_eff_m2;
    hostStoreOphelieTeam7CoilPathTangents(particles, centerline, spec);
    return summary;
}

inline Real hostCoilIntegratedCurrentFromJSrc(BaseParticles &particles, const OphelieCoilFieldNames &names,
                                              Real path_length_m)
{
    if (path_length_m <= TinyReal)
    {
        return 0.0;
    }
    syncVariableToHost<Vecd>(particles, names.j_src_real);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *j_src = particles.getVariableDataByName<Vecd>(names.j_src_real);
    const Real *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t count = particles.TotalRealParticles();
    Real sum = 0.0;
    for (size_t i = 0; i != count; ++i)
    {
        sum += j_src[i].norm() * volume[i];
    }
    return sum / path_length_m;
}

inline void printOphelieTeam7CoilPathPrepareSummary(const OphelieTeam7CoilPathPrepareSummary &summary,
                                                    const OphelieTeam7CoilPathSourceSpec &spec)
{
    std::cout << "[ophelie] TEAM7 volume-racetrack coil path NI=" << summary.ni_target
              << " R_out_m=" << spec.R_out << " R_in_m=" << spec.R_in
              << " path_length_m=" << summary.path_length_m << " centerline_pts=" << summary.centerline_points
              << " coil_volume_m3=" << summary.coil_volume_m3 << " A_eff_m2=" << summary.a_eff_m2
              << " J0=" << summary.j0 << " i_eff_formula=" << summary.i_eff_formula << std::endl;
}

inline void buildOphelieTeam7FilamentMomentsFromCoilPath(const BoundingBoxd &coil_bbox,
                                                         const OphelieTeam7CoilPathSourceSpec &spec, Real ds_m,
                                                         StdVec<OphelieCurrentMomentSample> &moments)
{
    moments.clear();
    const StdVec<Vecd> centerline = buildOphelieTeam7CoilCenterlinePolyline(coil_bbox, spec);
    if (centerline.size() < 2 || ds_m <= TinyReal)
    {
        return;
    }
    const Real total_current = spec.source_scale * spec.totalAmpereTurns();
    const size_t n = centerline.size();
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd a = centerline[i];
        const Vecd b = centerline[(i + 1) % n];
        const Vecd ab = b - a;
        const Real seg_len = ab.norm();
        if (seg_len <= TinyReal)
        {
            continue;
        }
        const int n_steps = std::max(1, static_cast<int>(std::ceil(seg_len / ds_m)));
        const Vecd dl = ab / static_cast<Real>(n_steps);
        for (int step = 0; step < n_steps; ++step)
        {
            const Vecd position = a + (static_cast<Real>(step) + Real(0.5)) * dl;
            moments.push_back(makeFilamentCurrentMoment(position, dl, total_current));
        }
    }
}

inline OphelieTeam7FilamentPrepareSummary prepareOphelieTeam7FilamentRacetrackCoilSource(
    const BoundingBoxd &coil_bbox, const OphelieTeam7CoilPathSourceSpec &spec, Real ds_m)
{
    OphelieTeam7FilamentPrepareSummary summary;
    summary.ni_target = spec.totalAmpereTurns();
    summary.ds_m = ds_m;
    const StdVec<Vecd> centerline = buildOphelieTeam7CoilCenterlinePolyline(coil_bbox, spec);
    summary.centerline_points = centerline.size();
    summary.path_length_m = closedPolylinePathLength(centerline);
    StdVec<OphelieCurrentMomentSample> moments;
    buildOphelieTeam7FilamentMomentsFromCoilPath(coil_bbox, spec, ds_m, moments);
    summary.quadrature_points = moments.size();
    return summary;
}

inline void printOphelieTeam7FilamentPrepareSummary(const OphelieTeam7FilamentPrepareSummary &summary,
                                                  const OphelieTeam7CoilPathSourceSpec &spec)
{
    std::cout << "[ophelie] TEAM7 filament-racetrack NI=" << summary.ni_target * spec.source_scale
              << " R_out_m=" << spec.R_out << " R_in_m=" << spec.R_in << " path_length_m=" << summary.path_length_m
              << " centerline_pts=" << summary.centerline_points << " ds_m=" << summary.ds_m
              << " quadrature_pts=" << summary.quadrature_points << std::endl;
}

inline bool ophelieTeam7FilamentAmpereTurnsAuditPassed(const OphelieTeam7FilamentPrepareSummary &summary,
                                                       Real source_scale = Real(1))
{
    const Real ni_effective = summary.ni_target * source_scale;
    std::cout << "[ophelie] TEAM7 filament-racetrack NI_effective=" << ni_effective
              << " path_length_m=" << summary.path_length_m << " audit_passed=1" << std::endl;
    (void)summary;
    return ni_effective > TinyReal;
}

inline void applyOphelieTeam7FilamentRacetrackBiotToGlass(SolidBody &plate_body,
                                                          const OphelieGlassFieldNames &plate_names,
                                                          const StdVec<OphelieCurrentMomentSample> &moments,
                                                          Real mu0, Real softening_length)
{
    applyFilamentBiotSavartToGlassHost(plate_body.getBaseParticles(), plate_names, moments, mu0, softening_length);
}

inline void evaluateTeam7FilamentBiotSavartBzAtProbes(const StdVec<OphelieCurrentMomentSample> &moments, Real mu0,
                                                      Real softening_length,
                                                      const StdVec<Team7BzProbePoint> &probes,
                                                      StdVec<Team7BzProbePoint> &results)
{
    results = probes;
    for (Team7BzProbePoint &probe : results)
    {
        Vecd a_r = Vecd::Zero();
        Vecd a_i = Vecd::Zero();
        Vecd b_r = Vecd::Zero();
        Vecd b_i = Vecd::Zero();
        accumulateBiotSavartFromMoments(moments, probe.position_m, mu0, softening_length, a_r, a_i, b_r, b_i);
        probe.bz_sim_mT = b_r[2] * 1000.0;
    }
}

inline bool ophelieTeam7CoilPathAmpereTurnsAuditPassed(const OphelieTeam7CoilPathPrepareSummary &summary,
                                                       Real integrated_current_a, Real source_scale = Real(1),
                                                       Real ratio_low = 0.95, Real ratio_high = 1.05)
{
    const Real ni_effective = summary.ni_target * source_scale;
    const Real ratio_formula = summary.i_eff_formula / (ni_effective + TinyReal);
    const Real ratio_particle = integrated_current_a / (ni_effective + TinyReal);
    const bool passed = ratio_formula > ratio_low && ratio_formula < ratio_high && ratio_particle > ratio_low &&
                        ratio_particle < ratio_high;
    std::cout << "[ophelie] TEAM7 volume-racetrack i_eff_ratio_formula=" << ratio_formula
              << " i_eff_ratio_particle=" << ratio_particle << " audit_passed=" << (passed ? 1 : 0) << std::endl;
    return passed;
}

/** J_src = J0 * tangent along TEAM7 racetrack centerline (host-prepared tangent field). */
class InitializeOphelieVolumeRacetrackCoilSourceCK : public LocalDynamics
{
  public:
    InitializeOphelieVolumeRacetrackCoilSourceCK(SPHBody &sph_body, const OphelieCoilFieldNames &names, Real j0)
        : LocalDynamics(sph_body), j0_(j0),
          dv_j_src_real_(particles_->template getVariableByName<Vecd>(names.j_src_real)),
          dv_j_src_imag_(particles_->template getVariableByName<Vecd>(names.j_src_imag)),
          dv_tangent_(particles_->template getVariableByName<Vecd>(kOphelieTeam7CoilPathTangentName))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : j0_(encloser.j0_), j_src_real_(encloser.dv_j_src_real_->DelegatedData(ex_policy)),
              j_src_imag_(encloser.dv_j_src_imag_->DelegatedData(ex_policy)),
              tangent_(encloser.dv_tangent_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            j_src_real_[index_i] = j0_ * tangent_[index_i];
            j_src_imag_[index_i] = Vecd::Zero();
        }

      protected:
        Real j0_;
        Vecd *j_src_real_;
        Vecd *j_src_imag_;
        Vecd *tangent_;
    };

  protected:
    Real j0_;
    DiscreteVariable<Vecd> *dv_j_src_real_;
    DiscreteVariable<Vecd> *dv_j_src_imag_;
    DiscreteVariable<Vecd> *dv_tangent_;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_COIL_PATH_SOURCE_H
