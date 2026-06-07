#ifndef APHI_TEAM7_NATIVE_COIL_SOURCE_HELPERS_H
#define APHI_TEAM7_NATIVE_COIL_SOURCE_HELPERS_H

#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"
#include "sphinxsys.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Official TEAM7 diagram radii [m]; STL bbox is used only for placement, not to infer R_out/R_in. */
struct Team7CoilPathSourceSpec
{
    Real turns = 2742.0;
    Real current_per_turn = 1.0;
    Real R_out = 0.050;
    Real R_in = 0.025;
    bool reverse_winding = false;
    bool use_particle_volume_cross_section = true;
    Real fallback_cross_section_area = 0.025 * 0.100;
    Real source_scale = 1.0;
    int polyline_segments_per_section = 32;

    Real totalAmpereTurns() const { return turns * current_per_turn; }
    Real centerlineRadius() const { return 0.5 * (R_out + R_in); }
};

inline constexpr const char *kTeam7CoilSourceTangentName = "Team7CoilSourceTangent";

struct Team7CoilBoundingBox
{
    Vecd lower = Vecd::Constant(std::numeric_limits<Real>::max());
    Vecd upper = Vecd::Constant(std::numeric_limits<Real>::lowest());
    bool valid = false;
};

inline StdVec<Vecd> buildTeam7CoilCenterlinePolyline(const Team7CoilBoundingBox &coil_bbox,
                                                     const Team7CoilPathSourceSpec &spec)
{
    StdVec<Vecd> pts;
    if (!coil_bbox.valid)
    {
        return pts;
    }

    const Real R_out = spec.R_out;
    const Real R_in = spec.R_in;
    const Real R_c = spec.centerlineRadius();
    const Real x_min = coil_bbox.lower[0];
    const Real x_max = coil_bbox.upper[0];
    const Real y_min = coil_bbox.lower[1];
    const Real y_max = coil_bbox.upper[1];
    const Real z_c = 0.5 * (coil_bbox.lower[2] + coil_bbox.upper[2]);

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

inline void hostStoreTeam7CoilSourceTangents(SPHBody &coil_body, const StdVec<Vecd> &centerline,
                                             const Team7CoilPathSourceSpec &spec)
{
    BaseParticles &particles = coil_body.getBaseParticles();
    particles.registerStateVariable<Vecd>(kTeam7CoilSourceTangentName, ZeroData<Vecd>::value);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, kTeam7CoilSourceTangentName);
    const size_t count = particles.TotalRealParticles();
    Vecd *tangents = particles.getVariableDataByName<Vecd>(kTeam7CoilSourceTangentName);
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
    syncVariableToDevice<Vecd>(particles, kTeam7CoilSourceTangentName);
}

/** rhs_a = J0 * tangent [A/m^2]; J0 = source_scale * NI / A_eff. */
class AssignTeam7CoilPathImpressedCurrentRhsCK : public LocalDynamics
{
  public:
    AssignTeam7CoilPathImpressedCurrentRhsCK(SPHBody &sph_body, const AphiBlockNames &rhs_block, Real current_density_j0)
        : LocalDynamics(sph_body),
          dv_rhs_a_real_(particles_->template getVariableByName<Vecd>(rhs_block.a_real)),
          dv_rhs_a_imag_(particles_->template getVariableByName<Vecd>(rhs_block.a_imag)),
          dv_tangent_(particles_->template getVariableByName<Vecd>(kTeam7CoilSourceTangentName)),
          current_density_j0_(current_density_j0)
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : rhs_a_real_(encloser.dv_rhs_a_real_->DelegatedData(ex_policy)),
              rhs_a_imag_(encloser.dv_rhs_a_imag_->DelegatedData(ex_policy)),
              tangent_(encloser.dv_tangent_->DelegatedData(ex_policy)), current_density_j0_(encloser.current_density_j0_)
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            rhs_a_real_[index_i] = current_density_j0_ * tangent_[index_i];
            rhs_a_imag_[index_i] = Vecd::Zero();
        }

      protected:
        Vecd *rhs_a_real_;
        Vecd *rhs_a_imag_;
        Vecd *tangent_;
        Real current_density_j0_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_rhs_a_real_;
    DiscreteVariable<Vecd> *dv_rhs_a_imag_;
    DiscreteVariable<Vecd> *dv_tangent_;
    Real current_density_j0_;
};

inline Real hostCoilVolume(SPHBody &coil_body)
{
    BaseParticles &particles = coil_body.getBaseParticles();
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

inline Real hostCoilIntegratedCurrentFromRhs(SPHBody &coil_body, const AphiBlockNames &rhs_block, Real path_length)
{
    if (path_length <= TinyReal)
    {
        return 0.0;
    }
    BaseParticles &particles = coil_body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, rhs_block.a_real);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *rhs = particles.getVariableDataByName<Vecd>(rhs_block.a_real);
    const Real *volume = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t count = particles.TotalRealParticles();
    Real sum = 0.0;
    for (size_t i = 0; i != count; ++i)
    {
        sum += rhs[i].norm() * volume[i];
    }
    return sum / path_length;
}

inline Real team7CoilPathCurrentDensityJ0(const Team7CoilBoundingBox &coil_bbox, const Team7CoilPathSourceSpec &spec,
                                          Real coil_volume)
{
    const StdVec<Vecd> centerline = buildTeam7CoilCenterlinePolyline(coil_bbox, spec);
    const Real path_length = closedPolylinePathLength(centerline);
    const Real a_eff = spec.use_particle_volume_cross_section && path_length > TinyReal
                           ? coil_volume / path_length
                           : spec.fallback_cross_section_area;
    return spec.source_scale * spec.totalAmpereTurns() / (a_eff + TinyReal);
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_NATIVE_COIL_SOURCE_HELPERS_H
