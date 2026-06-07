#ifndef APHI_TEAM7_NATIVE_MR_PROBE_HELPERS_H
#define APHI_TEAM7_NATIVE_MR_PROBE_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_geometry_config.h"
#include "electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline Real team7NativeProbeQuerySpacingM(const AphiTeam7NativeGeometryConfig &config)
{
    return config.dp_reference;
}

/** Probe h_ratio at finest resolution (framework: h_ref / h_local). */
inline Real team7NativeProbeHRatio(const AphiTeam7NativeGeometryConfig &config)
{
    return team7NativeFinestSmoothingLengthRatio(config);
}

/** Same support test as NeighborBuilderInnerAdaptive (framework neighborhood.cpp). */
inline bool team7NativeProbeWithinAdaptiveKernelSupport(const Kernel &kernel, Real distance, Real h_ratio_probe,
                                                        Real h_ratio_particle)
{
    const Real h_ratio_min = SMIN(h_ratio_probe, h_ratio_particle);
    return distance < kernel.CutOffRadius(h_ratio_min);
}

inline AphiProbeLineSample sampleNearestAirBWithTeam7KernelSupport(
    BaseParticles &air_particles, const Kernel &kernel, const AphiTeam7NativeGeometryConfig &geometry_config,
    const Vecd &sample_position, const char *b_real_name, const char *b_imag_name)
{
    AphiProbeLineSample sample;
    sample.position = sample_position;
    const Real h_ratio_probe = team7NativeProbeHRatio(geometry_config);
    syncVariableToHost<Real>(air_particles, "SmoothingLengthRatio");
    syncVariableToHost<Vecd>(air_particles, "Position");
    syncVariableToHost<Vecd>(air_particles, b_real_name);
    syncVariableToHost<Vecd>(air_particles, b_imag_name);
    const Vecd *positions = air_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *b_real = air_particles.getVariableDataByName<Vecd>(b_real_name);
    const Vecd *b_imag = air_particles.getVariableDataByName<Vecd>(b_imag_name);
    const Real *h_ratio = air_particles.getVariableDataByName<Real>("SmoothingLengthRatio");
    Real best_distance_squared = std::numeric_limits<Real>::max();
    const size_t total_real_particles = air_particles.TotalRealParticles();
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        const Real distance_squared = (positions[i] - sample_position).squaredNorm();
        const Real distance = std::sqrt(distance_squared);
        const Real h_ratio_particle = h_ratio[i];
        if (!team7NativeProbeWithinAdaptiveKernelSupport(kernel, distance, h_ratio_probe, h_ratio_particle))
        {
            continue;
        }
        if (distance_squared < best_distance_squared)
        {
            best_distance_squared = distance_squared;
            sample.nearest_particle_index = i;
            sample.neighbor_distance = distance;
            sample.b_real = b_real[i];
            sample.b_imag = b_imag[i];
            sample.b_abs = vec3Abs(sample.b_real);
        }
    }
    return sample;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_NATIVE_MR_PROBE_HELPERS_H
