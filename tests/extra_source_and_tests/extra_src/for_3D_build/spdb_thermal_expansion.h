/**
 * @file spdb_thermal_expansion.h
 * @brief Isotropic thermal stretch helper for programmed-temperature demo.
 * @details demo_assumption α. Scales current positions from a fixed origin and
 *          multiplies F by the same factor so free expansion stays consistent.
 *          Do not also apply a second thermal strain in the constitutive stress.
 */
#ifndef SPDB_THERMAL_EXPANSION_H
#define SPDB_THERMAL_EXPANSION_H

#include "sphinxsys.h"

namespace spdb
{

inline SPH::Real thermal_stretch_from_C(SPH::Real temperature_C, SPH::Real alpha_per_K,
                                        SPH::Real t_ref_C = 20.0)
{
    return 1.0 + alpha_per_K * (temperature_C - t_ref_C);
}

inline SPH::Real thermal_stretch_from_K(SPH::Real temperature_K, SPH::Real alpha_per_K,
                                        SPH::Real t_ref_K = 293.15)
{
    return 1.0 + alpha_per_K * (temperature_K - t_ref_K);
}

/** Scale positions about origin and F *= factor (isotropic). */
inline void apply_isotropic_thermal_factor(SPH::BaseParticles &particles, SPH::Real factor,
                                          const SPH::Vecd &origin)
{
    using SPH::Matd;
    using SPH::Real;
    using SPH::Vecd;
    if (!(factor > 0.0) || std::abs(factor - 1.0) < 1.0e-15)
        return;
    const size_t n = particles.TotalRealParticles();
    Vecd *pos = particles.ParticlePositions();
    Matd *F = particles.getVariableDataByName<Matd>("DeformationGradient");
    for (size_t i = 0; i < n; ++i)
    {
        pos[i] = origin + factor * (pos[i] - origin);
        F[i] = factor * F[i];
    }
}

inline void assign_uniform_temperature(SPH::BaseParticles &particles, SPH::Real temperature_K)
{
    SPH::Real *T = particles.getVariableDataByName<SPH::Real>("Temperature");
    for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
        T[i] = temperature_K;
}

inline void temperature_minmax(SPH::BaseParticles &particles, SPH::Real *tmin, SPH::Real *tmax)
{
    SPH::Real *T = particles.getVariableDataByName<SPH::Real>("Temperature");
    *tmin = T[0];
    *tmax = T[0];
    for (size_t i = 1; i < particles.TotalRealParticles(); ++i)
    {
        *tmin = std::min(*tmin, T[i]);
        *tmax = std::max(*tmax, T[i]);
    }
}

} // namespace spdb

#endif
