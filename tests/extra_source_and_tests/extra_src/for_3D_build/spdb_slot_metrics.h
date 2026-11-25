/**
 * @file spdb_slot_metrics.h
 * @brief Host-side cavity / join-face metrics for the SPDB blind-slot demo.
 * @details No claimed cavity volume. Net height and widths are particle-center
 *          reconstructions. Contact penetration must stay separate from creep.
 */
#ifndef SPDB_SLOT_METRICS_H
#define SPDB_SLOT_METRICS_H

#include "spdb_geometry.h"

#include <cmath>
#include <limits>
#include <vector>

namespace SPH
{
// Keep helpers in the global `spdb` namespace used by Arrhenius (not SPH::spdb).
} // namespace SPH

namespace spdb
{

struct SlotCavityMetrics
{
    size_t floor_n = 0;
    size_t roof_n = 0;
    size_t contact_lower_n = 0;
    size_t contact_upper_n = 0;
    SPH::Real center_net_height = 0.0;   // mean roof z - mean floor z
    SPH::Real min_clearance = 0.0;       // min roof-floor pair with matching xy
    SPH::Real width_x = 0.0;             // floor particle span in x
    SPH::Real width_y = 0.0;
    SPH::Real roof_mean_z = 0.0;
    SPH::Real floor_mean_z = 0.0;
    SPH::Real contact_ring_upper_mean_z = 0.0;
    SPH::Real top_wall_deflection = 0.0; // roof_mean_z - contact_ring_upper_mean_z
};

inline bool in_slot_xy(const SPH::Vecd &p, SPH::Real margin = 0.0)
{
    using SPH::slot_half_xy;
    return std::abs(p[0]) < slot_half_xy - margin && std::abs(p[1]) < slot_half_xy - margin;
}

inline bool outside_slot_xy(const SPH::Vecd &p, SPH::Real margin = 0.0)
{
    using SPH::slot_half_xy;
    return std::abs(p[0]) >= slot_half_xy + margin || std::abs(p[1]) >= slot_half_xy + margin;
}

inline SPH::IndexVector ids_near_z_in_slot(SPH::Vecd *pos, size_t n, SPH::Real z_target, SPH::Real z_tol)
{
    using SPH::particle_spacing_ref;
    SPH::IndexVector ids;
    for (size_t i = 0; i < n; ++i)
        if (std::abs(pos[i][2] - z_target) <= z_tol && in_slot_xy(pos[i], -0.25 * particle_spacing_ref))
            ids.push_back(i);
    return ids;
}

inline SPH::IndexVector ids_near_z_outside_slot(SPH::Vecd *pos, size_t n, SPH::Real z_target, SPH::Real z_tol)
{
    SPH::IndexVector ids;
    for (size_t i = 0; i < n; ++i)
        if (std::abs(pos[i][2] - z_target) <= z_tol && outside_slot_xy(pos[i], 0.0))
            ids.push_back(i);
    return ids;
}

inline SlotCavityMetrics measure_slot_cavity(SPH::Vecd *pos_l, size_t n_l, SPH::Vecd *pos_u, size_t n_u, SPH::Real dp)
{
    using SPH::Real;
    using SPH::IndexVector;
    using SPH::slot_depth;
    SlotCavityMetrics m;
    // Blind slot removes z in [-slot_depth, 0]. Floor particles sit half a spacing
    // inside the remaining solid, so their centers are near -slot_depth - 0.5*dp.
    // Upper roof over the cavity sits at +0.5*dp inside the upper body.
    const Real z_floor = -slot_depth - 0.5 * dp;
    const Real z_join_l = -0.5 * dp;
    const Real z_join_u = 0.5 * dp;
    IndexVector floor = ids_near_z_in_slot(pos_l, n_l, z_floor, 0.55 * dp);
    IndexVector roof = ids_near_z_in_slot(pos_u, n_u, z_join_u, 0.55 * dp);
    IndexVector contact_l = ids_near_z_outside_slot(pos_l, n_l, z_join_l, 0.55 * dp);
    IndexVector contact_u = ids_near_z_outside_slot(pos_u, n_u, z_join_u, 0.55 * dp);
    m.floor_n = floor.size();
    m.roof_n = roof.size();
    m.contact_lower_n = contact_l.size();
    m.contact_upper_n = contact_u.size();
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
            if (std::sqrt(dx * dx + dy * dy) > 0.75 * dp)
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

#endif // SPDB_SLOT_METRICS_H
