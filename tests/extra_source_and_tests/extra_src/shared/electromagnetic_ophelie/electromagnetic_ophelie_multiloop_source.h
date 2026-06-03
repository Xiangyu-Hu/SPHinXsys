#ifndef ELECTROMAGNETIC_OPHELIE_MULTILOOP_SOURCE_H
#define ELECTROMAGNETIC_OPHELIE_MULTILOOP_SOURCE_H

#include "electromagnetic_ophelie_current_moment_source.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieCircularLoopSpec
{
    Vecd center = Vecd::Zero();
    Real radius = 0.1;
    Real current = 1.0;
    size_t segments = 360;
};

struct OphelieMultiloopCoilSpec
{
    Vecd stack_center = Vecd::Zero();
    Real loop_radius = 0.1;
    Real z_min = 0.0;
    Real z_max = 0.0;
    size_t num_loops = 1;
    Real current_per_loop = 1.0;
    size_t segments_per_loop = 360;
};

inline Real analyticCircularLoopBzOnAxis(Real z_observer, Real z_loop, Real radius, Real current, Real mu0)
{
    const Real dz = z_observer - z_loop;
    const Real denom = std::pow(radius * radius + dz * dz, 1.5);
    return mu0 * current * radius * radius / (2.0 * denom);
}

inline Real analyticMultiloopBzOnAxis(Real z_observer, const OphelieMultiloopCoilSpec &spec, Real mu0)
{
    if (spec.num_loops == 0)
    {
        return 0.0;
    }
    Real bz = 0.0;
    const Real dz_stack = spec.num_loops > 1 ? (spec.z_max - spec.z_min) / static_cast<Real>(spec.num_loops - 1) : 0.0;
    for (size_t k = 0; k < spec.num_loops; ++k)
    {
        const Real z_loop = spec.num_loops > 1 ? spec.z_min + static_cast<Real>(k) * dz_stack : 0.5 * (spec.z_min + spec.z_max);
        bz += analyticCircularLoopBzOnAxis(z_observer, z_loop, spec.loop_radius, spec.current_per_loop, mu0);
    }
    return bz;
}

inline void buildCircularLoopFilamentMoments(const OphelieCircularLoopSpec &loop,
                                             StdVec<OphelieCurrentMomentSample> &moments)
{
    moments.clear();
    moments.reserve(loop.segments);
    const Real dtheta = 2.0 * Pi / static_cast<Real>(loop.segments);
    for (size_t i = 0; i < loop.segments; ++i)
    {
        const Real theta = static_cast<Real>(i) * dtheta;
        const Real theta_next = static_cast<Real>(i + 1) * dtheta;
        const Vecd p0(loop.center[0] + loop.radius * std::cos(theta), loop.center[1] + loop.radius * std::sin(theta),
                      loop.center[2]);
        const Vecd p1(loop.center[0] + loop.radius * std::cos(theta_next),
                      loop.center[1] + loop.radius * std::sin(theta_next), loop.center[2]);
        const Vecd midpoint = 0.5 * (p0 + p1);
        const Vecd dl = p1 - p0;
        moments.push_back(makeFilamentCurrentMoment(midpoint, dl, loop.current));
    }
}

inline void buildMultiloopFilamentMoments(const OphelieMultiloopCoilSpec &spec,
                                          StdVec<OphelieCurrentMomentSample> &moments)
{
    moments.clear();
    if (spec.num_loops == 0)
    {
        return;
    }
    const Real dz_stack =
        spec.num_loops > 1 ? (spec.z_max - spec.z_min) / static_cast<Real>(spec.num_loops - 1) : 0.0;
    for (size_t k = 0; k < spec.num_loops; ++k)
    {
        const Real z_loop =
            spec.num_loops > 1 ? spec.z_min + static_cast<Real>(k) * dz_stack : 0.5 * (spec.z_min + spec.z_max);
        OphelieCircularLoopSpec loop;
        loop.center = Vecd(spec.stack_center[0], spec.stack_center[1], z_loop);
        loop.radius = spec.loop_radius;
        loop.current = spec.current_per_loop;
        loop.segments = spec.segments_per_loop;
        StdVec<OphelieCurrentMomentSample> loop_moments;
        buildCircularLoopFilamentMoments(loop, loop_moments);
        moments.insert(moments.end(), loop_moments.begin(), loop_moments.end());
    }
}

inline Real filamentBzOnAxisFromMoments(const StdVec<OphelieCurrentMomentSample> &moments, const Vecd &observer,
                                        Real mu0, Real softening_length = 0.0)
{
    Vecd a_r = Vecd::Zero();
    Vecd a_i = Vecd::Zero();
    Vecd b_r = Vecd::Zero();
    Vecd b_i = Vecd::Zero();
    accumulateBiotSavartFromMoments(moments, observer, mu0, softening_length, a_r, a_i, b_r, b_i);
    return b_r[2];
}

inline void applyMultiloopFilamentBiotToGlass(BaseParticles &glass_particles, const OphelieGlassFieldNames &glass_names,
                                              const OphelieMultiloopCoilSpec &spec, Real mu0, Real softening_length)
{
    StdVec<OphelieCurrentMomentSample> moments;
    buildMultiloopFilamentMoments(spec, moments);
    applyFilamentBiotSavartToGlassHost(glass_particles, glass_names, moments, mu0, softening_length);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_MULTILOOP_SOURCE_H
