#ifndef APHI_EM_OBSERVABLE_HELPERS_H
#define APHI_EM_OBSERVABLE_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <cmath>
#include <functional>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline bool isInsideTeam7PhysicalBox(const Vecd &position, Real body_length, Real body_height, Real body_width)
{
    return position[0] >= 0.0 && position[0] <= body_length && position[1] >= 0.0 && position[1] <= body_height &&
           position[2] >= 0.0 && position[2] <= body_width;
}

inline bool isInsidePassiveAirShellParticle(const Vecd &position, Real body_length, Real body_height, Real body_width)
{
    return !isInsideTeam7PhysicalBox(position, body_length, body_height, body_width);
}

inline bool team7ParticleInPhysicalRegion(const Vecd &position, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                          Real body_length, Real body_height, Real body_width,
                                          AphiBenchmarkMaterialRegion region)
{
    return isInsideTeam7PhysicalBox(position, body_length, body_height, body_width) &&
           team7ParticleInRegion(position, layout, region);
}

struct AphiRegionalElectromagneticObservables
{
    Real E_real_max = 0.0;
    Real E_imag_max = 0.0;
    Real E_L2 = 0.0;
    Real J_real_max = 0.0;
    Real J_imag_max = 0.0;
    Real J_L2 = 0.0;
    Real Joule_power = 0.0;
    Real Joule_max = 0.0;
    Real Joule_L2 = 0.0;
    Real solution_block_max = 0.0;
};

inline Real hostParticleRegionVolWeightedVecdL2(
    BaseParticles &particles, const std::string &variable_name, const Vecd *positions, size_t total_real_particles,
    const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Vecd>(particles, variable_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *values = particles.getVariableDataByName<Vecd>(variable_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        sum_squared += vol[i] * values[i].squaredNorm();
    }
    return std::sqrt(sum_squared);
}

inline Real hostParticleRegionVolWeightedScalarL2(
    BaseParticles &particles, const std::string &variable_name, const Vecd *positions, size_t total_real_particles,
    const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        sum_squared += vol[i] * values[i] * values[i];
    }
    return std::sqrt(sum_squared);
}

inline Real hostParticleRegionVecdMax(BaseParticles &particles, const std::string &variable_name, const Vecd *positions,
                                    size_t total_real_particles,
                                    const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Vecd>(particles, variable_name);
    const Vecd *values = particles.getVariableDataByName<Vecd>(variable_name);
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        max_value = std::max(max_value, values[i].norm());
    }
    return max_value;
}

inline Real hostParticleRegionScalarMax(BaseParticles &particles, const std::string &variable_name,
                                        const Vecd *positions, size_t total_real_particles,
                                        const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Real>(particles, variable_name);
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        max_value = std::max(max_value, values[i]);
    }
    return max_value;
}

inline Real hostParticleRegionVolWeightedJoulePower(BaseParticles &particles, const Vecd *positions,
                                                     size_t total_real_particles,
                                                     const std::function<bool(const Vecd &)> &in_region,
                                                     const std::string &joule_source_name)
{
    return hostVolWeightedSum(
        particles, joule_source_name, total_real_particles,
        [&](const Vecd &position) { return in_region(position); }, positions);
}

/** E/J/Joule observables from post-processed Joule fields on a particle region. */
inline AphiRegionalElectromagneticObservables hostRegionalElectromagneticObservablesFromJouleFields(
    BaseParticles &particles, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_names,
    const Vecd *positions, size_t total_real_particles, Real body_length, Real body_height, Real body_width,
    Real core_shell, const std::function<bool(const Vecd &)> &in_region, Real sigma_for_current_density)
{
    AphiRegionalElectromagneticObservables observables;
    observables.E_real_max =
        hostParticleRegionVecdMax(particles, joule_names.electric_field_a_real, positions, total_real_particles, in_region);
    observables.E_imag_max =
        hostParticleRegionVecdMax(particles, joule_names.electric_field_a_imag, positions, total_real_particles, in_region);
    observables.E_L2 = hostParticleRegionVolWeightedVecdL2(particles, joule_names.electric_field_a_real, positions,
                                                           total_real_particles, in_region);
    observables.J_real_max = sigma_for_current_density * observables.E_real_max;
    observables.J_imag_max = sigma_for_current_density * observables.E_imag_max;
    observables.J_L2 = sigma_for_current_density * observables.E_L2;
    observables.Joule_power = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles,
                                                                      in_region, joule_names.joule_heat_source);
    observables.Joule_max = hostParticleRegionScalarMax(particles, joule_names.joule_heat_source, positions,
                                                        total_real_particles, in_region);
    observables.Joule_L2 = hostParticleRegionVolWeightedScalarL2(particles, joule_names.joule_heat_source, positions,
                                                                 total_real_particles, in_region);
    observables.solution_block_max =
        coreBlockMaxAbsNormInRegionPredicate(particles, names.solution, positions, total_real_particles, body_length,
                                             body_height, body_width, core_shell, in_region);
    return observables;
}

inline AphiRegionalElectromagneticObservables hostRegionalElectromagneticObservablesFromTeam7Region(
    BaseParticles &particles, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_names,
    const Vecd *positions, size_t total_real_particles, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
    AphiBenchmarkMaterialRegion region, Real body_length, Real body_height, Real body_width, Real core_shell)
{
    const auto in_region = [&](const Vecd &position) {
        return team7ParticleInRegion(position, layout, region);
    };
    Real sigma = layout.air.sigma;
    if (region == AphiBenchmarkMaterialRegion::Conductor)
    {
        sigma = layout.conductor_material.sigma;
    }
    else if (region == AphiBenchmarkMaterialRegion::Coil)
    {
        sigma = layout.coil_material.sigma;
    }
    return hostRegionalElectromagneticObservablesFromJouleFields(particles, names, joule_names, positions,
                                                                 total_real_particles, body_length, body_height,
                                                                 body_width, core_shell, in_region, sigma);
}

struct AphiConductorInterfaceSpikeMetrics
{
    size_t interface_particle_count = 0;
    size_t bulk_conductor_particle_count = 0;
    Real interface_J_max = 0.0;
    Real mean_J_interface = 0.0;
    Real bulk_J_max = 0.0;
    Real mean_J_bulk = 0.0;
    Real interface_Joule_max = 0.0;
    Real mean_Joule_interface = 0.0;
    Real bulk_Joule_max = 0.0;
    Real mean_Joule_bulk = 0.0;
    Real J_spike_ratio_max = 0.0;
    Real J_spike_ratio_mean = 0.0;
    Real Joule_spike_ratio_max = 0.0;
    Real Joule_spike_ratio_mean = 0.0;
    bool spike_warning = false;
    bool spike_hard_pass = false;
    bool has_bulk_reference = false;
};

inline Real hostRegionMeanScalar(BaseParticles &particles, const std::string &variable_name, const Vecd *positions,
                                 size_t total_real_particles, const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Real>(particles, variable_name);
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    Real sum = 0.0;
    size_t count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        sum += values[i];
        ++count;
    }
    return sum / (static_cast<Real>(count) + TinyReal);
}

inline Real hostRegionMeanVecdNorm(BaseParticles &particles, const std::string &variable_name, const Vecd *positions,
                                   size_t total_real_particles, const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Vecd>(particles, variable_name);
    const Vecd *values = particles.getVariableDataByName<Vecd>(variable_name);
    Real sum = 0.0;
    size_t count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        sum += values[i].norm();
        ++count;
    }
    return sum / (static_cast<Real>(count) + TinyReal);
}

/** Stage 10.14-B: interface-band vs bulk conductor spike metrics (warning @10, hard fail @50). */
inline AphiConductorInterfaceSpikeMetrics hostConductorInterfaceSpikeMetrics(
    BaseParticles &particles, const AphiJouleHeatingFieldNames &joule_names, const Vecd *positions,
    size_t total_real_particles, const benchmark::AphiTeam7LikeUnitBoxLayout &layout, Real body_length,
    Real body_height, Real body_width, Real core_shell, Real x_interface, Real interface_band_half_width, Real dp,
    Real warning_spike = 10.0, Real hard_spike = 50.0)
{
    const Real sigma = layout.conductor_material.sigma;
    const auto in_conductor = [&](const Vecd &position) {
        return isInsideTeam7PhysicalBox(position, body_length, body_height, body_width) &&
               team7ParticleInRegion(position, layout, AphiBenchmarkMaterialRegion::Conductor);
    };
    const auto interface_band = [&](const Vecd &position) {
        return in_conductor(position) && std::abs(position[0] - x_interface) <= interface_band_half_width;
    };
    const auto bulk_conductor = [&](const Vecd &position) {
        return in_conductor(position) &&
               std::abs(position[0] - x_interface) > interface_band_half_width + TinyReal;
    };

    AphiConductorInterfaceSpikeMetrics metrics;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (interface_band(positions[i]))
        {
            metrics.interface_particle_count += 1;
        }
        if (bulk_conductor(positions[i]))
        {
            metrics.bulk_conductor_particle_count += 1;
        }
    }

    metrics.interface_J_max =
        sigma * hostParticleRegionVecdMax(particles, joule_names.electric_field_a_real, positions, total_real_particles,
                                        interface_band);
    metrics.bulk_J_max =
        sigma * hostParticleRegionVecdMax(particles, joule_names.electric_field_a_real, positions, total_real_particles,
                                         bulk_conductor);
    metrics.mean_J_interface =
        sigma * hostRegionMeanVecdNorm(particles, joule_names.electric_field_a_real, positions, total_real_particles,
                                       interface_band);
    metrics.mean_J_bulk =
        sigma * hostRegionMeanVecdNorm(particles, joule_names.electric_field_a_real, positions, total_real_particles,
                                       bulk_conductor);
    metrics.interface_Joule_max = hostParticleRegionScalarMax(particles, joule_names.joule_heat_source, positions,
                                                            total_real_particles, interface_band);
    metrics.bulk_Joule_max = hostParticleRegionScalarMax(particles, joule_names.joule_heat_source, positions,
                                                       total_real_particles, bulk_conductor);
    metrics.mean_Joule_interface = hostRegionMeanScalar(particles, joule_names.joule_heat_source, positions,
                                                      total_real_particles, interface_band);
    metrics.mean_Joule_bulk = hostRegionMeanScalar(particles, joule_names.joule_heat_source, positions,
                                                   total_real_particles, bulk_conductor);

    metrics.has_bulk_reference =
        metrics.bulk_conductor_particle_count > 0 && metrics.bulk_J_max > TinyReal && metrics.bulk_Joule_max > TinyReal;
    if (!metrics.has_bulk_reference)
    {
        metrics.spike_warning = false;
        metrics.spike_hard_pass = false;
        return metrics;
    }

    metrics.J_spike_ratio_max = metrics.interface_J_max / (metrics.bulk_J_max + TinyReal);
    metrics.Joule_spike_ratio_max = metrics.interface_Joule_max / (metrics.bulk_Joule_max + TinyReal);
    metrics.J_spike_ratio_mean = metrics.mean_J_interface / (metrics.mean_J_bulk + TinyReal);
    metrics.Joule_spike_ratio_mean = metrics.mean_Joule_interface / (metrics.mean_Joule_bulk + TinyReal);
    metrics.spike_warning =
        metrics.J_spike_ratio_max > warning_spike || metrics.Joule_spike_ratio_max > warning_spike;
    metrics.spike_hard_pass =
        metrics.J_spike_ratio_max < hard_spike && metrics.Joule_spike_ratio_max < hard_spike;
    return metrics;
}

inline void printConductorInterfaceSpikeMetrics(const char *prefix, const AphiConductorInterfaceSpikeMetrics &metrics,
                                                Real warning_spike, Real hard_spike)
{
    std::cout << prefix << " interface_particles=" << metrics.interface_particle_count
              << " bulk_conductor_particles=" << metrics.bulk_conductor_particle_count
              << " interface_J_max=" << metrics.interface_J_max << " mean_J_interface=" << metrics.mean_J_interface
              << " bulk_J_max=" << metrics.bulk_J_max << " mean_J_bulk=" << metrics.mean_J_bulk
              << " interface_Joule_max=" << metrics.interface_Joule_max
              << " mean_Joule_interface=" << metrics.mean_Joule_interface
              << " bulk_Joule_max=" << metrics.bulk_Joule_max << " mean_Joule_bulk=" << metrics.mean_Joule_bulk
              << " J_spike_ratio_max=" << metrics.J_spike_ratio_max
              << " J_spike_ratio_mean=" << metrics.J_spike_ratio_mean
              << " Joule_spike_ratio_max=" << metrics.Joule_spike_ratio_max
              << " Joule_spike_ratio_mean=" << metrics.Joule_spike_ratio_mean
              << " has_bulk_reference=" << (metrics.has_bulk_reference ? 1 : 0)
              << " spike_hard_pass=" << (metrics.spike_hard_pass ? 1 : 0) << " warning_spike=" << warning_spike
              << " hard_spike=" << hard_spike;
    if (metrics.spike_warning)
    {
        std::cout << " warning=interface_spike";
    }
    if (!metrics.has_bulk_reference)
    {
        std::cout << " warning=no_bulk_reference_for_spike";
    }
    std::cout << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_EM_OBSERVABLE_HELPERS_H
