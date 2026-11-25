/**
 * @file spdb_dual_slot_metrics.h
 * @brief Particle-center cavity metrics for facing dual blind slots (nominal 20 mm height).
 * @details No claimed cavity volume. Penetration must stay separate from creep/thermal.
 */
#ifndef SPDB_DUAL_SLOT_METRICS_H
#define SPDB_DUAL_SLOT_METRICS_H

#include "sphinxsys.h"

#include <cmath>
#include <limits>
#include <vector>

namespace spdb
{

inline constexpr SPH::Real dual_slot_half_xy = 0.005; // 10 mm square
inline constexpr SPH::Real dual_slot_each_depth = 0.010;

struct DualSlotCavityMetrics
{
    size_t floor_n = 0;
    size_t roof_n = 0;
    size_t contact_lower_n = 0;
    size_t contact_upper_n = 0;
    SPH::Real center_net_height = 0.0;
    SPH::Real min_clearance = 0.0;
    SPH::Real width_x = 0.0;
    SPH::Real width_y = 0.0;
    SPH::Real floor_mean_z = 0.0;
    SPH::Real roof_mean_z = 0.0;
    SPH::Real contact_ring_upper_mean_z = 0.0;
    SPH::Real top_wall_deflection = 0.0;
    SPH::Real lower_height = 0.0;
    SPH::Real upper_height = 0.0;
    SPH::Real stack_height = 0.0;
};

inline bool dual_in_slot_xy(const SPH::Vecd &p, SPH::Real margin = 0.0)
{
    return std::abs(p[0]) < dual_slot_half_xy - margin && std::abs(p[1]) < dual_slot_half_xy - margin;
}

inline bool dual_outside_slot_xy(const SPH::Vecd &p, SPH::Real margin = 0.0)
{
    return std::abs(p[0]) >= dual_slot_half_xy + margin || std::abs(p[1]) >= dual_slot_half_xy + margin;
}

inline SPH::IndexVector dual_ids_near_z_in_slot(SPH::Vecd *pos, size_t n, SPH::Real z_target, SPH::Real z_tol,
                                               SPH::Real xy_margin)
{
    SPH::IndexVector ids;
    for (size_t i = 0; i < n; ++i)
        if (std::abs(pos[i][2] - z_target) <= z_tol && dual_in_slot_xy(pos[i], xy_margin))
            ids.push_back(i);
    return ids;
}

inline SPH::IndexVector dual_ids_near_z_outside_slot(SPH::Vecd *pos, size_t n, SPH::Real z_target, SPH::Real z_tol)
{
    SPH::IndexVector ids;
    for (size_t i = 0; i < n; ++i)
        if (std::abs(pos[i][2] - z_target) <= z_tol && dual_outside_slot_xy(pos[i], 0.0))
            ids.push_back(i);
    return ids;
}

inline void body_z_span(SPH::Vecd *pos, size_t n, SPH::Real *zmin, SPH::Real *zmax)
{
    *zmin = std::numeric_limits<SPH::Real>::max();
    *zmax = -std::numeric_limits<SPH::Real>::max();
    for (size_t i = 0; i < n; ++i)
    {
        *zmin = std::min(*zmin, pos[i][2]);
        *zmax = std::max(*zmax, pos[i][2]);
    }
}

/** In-slot particles within z_tol of the body extremum (floor=min z, roof=max z). */
inline SPH::IndexVector dual_ids_extremal_in_slot(SPH::Vecd *pos, size_t n, bool take_min_z, SPH::Real z_tol,
                                                 SPH::Real xy_margin)
{
    using SPH::Real;
    using SPH::IndexVector;
    IndexVector ids;
    Real z_ext = take_min_z ? std::numeric_limits<Real>::max() : -std::numeric_limits<Real>::max();
    bool any = false;
    for (size_t i = 0; i < n; ++i)
    {
        if (!dual_in_slot_xy(pos[i], xy_margin))
            continue;
        any = true;
        z_ext = take_min_z ? std::min(z_ext, pos[i][2]) : std::max(z_ext, pos[i][2]);
    }
    if (!any)
        return ids;
    for (size_t i = 0; i < n; ++i)
        if (dual_in_slot_xy(pos[i], xy_margin) && std::abs(pos[i][2] - z_ext) <= z_tol)
            ids.push_back(i);
    return ids;
}

/**
 * Dual facing slots: prefer nominal absolute z bands; if empty (bodies drifted),
 * fall back to in-slot z extrema so GPU/thermal runs still report a net height.
 */
inline DualSlotCavityMetrics measure_dual_slot_cavity(SPH::Vecd *pos_l, size_t n_l, SPH::Vecd *pos_u, size_t n_u,
                                                     SPH::Real dp)
{
    using SPH::Real;
    using SPH::IndexVector;
    DualSlotCavityMetrics m;
    const Real z_floor = -dual_slot_each_depth - 0.5 * dp;
    const Real z_roof = dual_slot_each_depth + 0.5 * dp;
    const Real z_join_l = -0.5 * dp;
    const Real z_join_u = 0.5 * dp;
    IndexVector floor = dual_ids_near_z_in_slot(pos_l, n_l, z_floor, 0.6 * dp, -0.25 * dp);
    IndexVector roof = dual_ids_near_z_in_slot(pos_u, n_u, z_roof, 0.6 * dp, -0.25 * dp);
    // Blind pockets: lower cavity floor = highest solid in slot XY; upper roof = lowest solid.
    // (min-z lower / max-z upper would hit outer bottoms/tops and inflate net height.)
    if (floor.empty())
        floor = dual_ids_extremal_in_slot(pos_l, n_l, /*take_min_z=*/false, 0.6 * dp, -0.25 * dp);
    if (roof.empty())
        roof = dual_ids_extremal_in_slot(pos_u, n_u, /*take_min_z=*/true, 0.6 * dp, -0.25 * dp);
    IndexVector contact_l = dual_ids_near_z_outside_slot(pos_l, n_l, z_join_l, 0.55 * dp);
    IndexVector contact_u = dual_ids_near_z_outside_slot(pos_u, n_u, z_join_u, 0.55 * dp);
    if (contact_l.empty() || contact_u.empty())
    {
        // Relative join-face bands from current body spans (outside slot XY).
        Real zl0 = 0.0, zl1 = 0.0, zu0 = 0.0, zu1 = 0.0;
        body_z_span(pos_l, n_l, &zl0, &zl1);
        body_z_span(pos_u, n_u, &zu0, &zu1);
        if (contact_l.empty())
            contact_l = dual_ids_near_z_outside_slot(pos_l, n_l, zl1 - 0.5 * dp, 0.75 * dp);
        if (contact_u.empty())
            contact_u = dual_ids_near_z_outside_slot(pos_u, n_u, zu0 + 0.5 * dp, 0.75 * dp);
    }
    m.floor_n = floor.size();
    m.roof_n = roof.size();
    m.contact_lower_n = contact_l.size();
    m.contact_upper_n = contact_u.size();

    Real zl0 = 0.0, zl1 = 0.0, zu0 = 0.0, zu1 = 0.0;
    body_z_span(pos_l, n_l, &zl0, &zl1);
    body_z_span(pos_u, n_u, &zu0, &zu1);
    m.lower_height = zl1 - zl0;
    m.upper_height = zu1 - zu0;
    m.stack_height = zu1 - zl0;

    if (floor.empty() || roof.empty())
        return m;

    Real floor_sum = 0.0;
    Real roof_sum = 0.0;
    Real xmin = std::numeric_limits<Real>::max();
    Real xmax = -std::numeric_limits<Real>::max();
    Real ymin = std::numeric_limits<Real>::max();
    Real ymax = -std::numeric_limits<Real>::max();
    for (size_t i : floor)
    {
        floor_sum += pos_l[i][2];
        xmin = std::min(xmin, pos_l[i][0]);
        xmax = std::max(xmax, pos_l[i][0]);
        ymin = std::min(ymin, pos_l[i][1]);
        ymax = std::max(ymax, pos_l[i][1]);
    }
    for (size_t j : roof)
        roof_sum += pos_u[j][2];
    m.floor_mean_z = floor_sum / static_cast<Real>(floor.size());
    m.roof_mean_z = roof_sum / static_cast<Real>(roof.size());
    m.center_net_height = m.roof_mean_z - m.floor_mean_z;
    m.width_x = xmax - xmin;
    m.width_y = ymax - ymin;

    Real min_clear = std::numeric_limits<Real>::max();
    for (size_t j : roof)
    {
        Real best = std::numeric_limits<Real>::max();
        for (size_t i : floor)
        {
            const Real dx = pos_u[j][0] - pos_l[i][0];
            const Real dy = pos_u[j][1] - pos_l[i][1];
            if (std::sqrt(dx * dx + dy * dy) > 0.85 * dp)
                continue;
            best = std::min(best, pos_u[j][2] - pos_l[i][2]);
        }
        if (best < std::numeric_limits<Real>::max() / 2.0)
            min_clear = std::min(min_clear, best);
    }
    m.min_clearance = (min_clear < std::numeric_limits<Real>::max() / 2.0) ? min_clear : m.center_net_height;

    if (!contact_u.empty())
    {
        Real sum = 0.0;
        for (size_t j : contact_u)
            sum += pos_u[j][2];
        m.contact_ring_upper_mean_z = sum / static_cast<Real>(contact_u.size());
        m.top_wall_deflection = m.roof_mean_z - m.contact_ring_upper_mean_z;
    }
    return m;
}

} // namespace spdb

#endif
