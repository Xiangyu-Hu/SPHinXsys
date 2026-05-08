#ifndef ELECTROMAGNETIC_MULTITURN_COIL_DRIVE_HPP
#define ELECTROMAGNETIC_MULTITURN_COIL_DRIVE_HPP

#include "electromagnetic_multiturn_coil_drive.h"

namespace SPH
{
namespace extra_electromagnetics
{
inline Vecd NormalizeVectorOrDefault(const Vecd &candidate, const Vecd &fallback)
{
    Real norm = candidate.norm();
    if (norm <= TinyReal)
    {
        return fallback;
    }
    return candidate / norm;
}

inline Vecd TangentialDirectionAroundAxis(const Vecd &position, const Vecd &center, const Vecd &axis)
{
    Vecd relative = position - center;
    Vecd radial = relative - relative.dot(axis) * axis;
    Real radial_norm = radial.norm();
    if (radial_norm <= TinyReal)
    {
        return ZeroData<Vecd>::value;
    }
    return axis.cross(radial) / (radial_norm + TinyReal);
}

inline MultiturnCoilDriveState BuildMultiturnCoilDriveState(BaseParticles &coil_particles,
                                                            const BoundingBoxd &coil_bounds,
                                                            const Vecd &coil_center,
                                                            Real angular_frequency,
                                                            const MultiturnCoilDriveConfig &config)
{
    MultiturnCoilDriveState state;
    state.source_direction_used =
        NormalizeVectorOrDefault(config.source_direction, Vecd(0.0, 0.0, 1.0));
    state.source_axis_used =
        NormalizeVectorOrDefault(config.source_axis, Vecd(0.0, 0.0, 1.0));
    state.harmonic_source_amplitude_used = config.fallback_harmonic_amplitude;
    state.frequency_source_real = config.fallback_harmonic_amplitude;
    state.frequency_source_imag = config.imag_scale * config.fallback_harmonic_amplitude;

    Vecd *positions = coil_particles.getVariableDataByName<Vecd>("Position");
    Real *volumetric_measure = coil_particles.getVariableDataByName<Real>("VolumetricMeasure");
    size_t total_real_particles = coil_particles.TotalRealParticles();
    Real coil_radius_volume_moment = 0.0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        state.volume_total_geom += volumetric_measure[i];
        if (config.use_circular_source)
        {
            Vecd radial = positions[i] - coil_center;
            radial -= radial.dot(state.source_axis_used) * state.source_axis_used;
            coil_radius_volume_moment += volumetric_measure[i] * radial.norm();
        }
    }

    if (config.use_circular_source)
    {
        state.mean_radius_geom =
            coil_radius_volume_moment / (state.volume_total_geom + TinyReal);
        state.projection_length_geom = 2.0 * Pi * state.mean_radius_geom;
    }
    else
    {
        Vecd coil_extent = coil_bounds.upper_ - coil_bounds.lower_;
        state.projection_length_geom =
            fabs(state.source_direction_used[0]) * coil_extent[0] +
            fabs(state.source_direction_used[1]) * coil_extent[1] +
            fabs(state.source_direction_used[2]) * coil_extent[2];
    }

    Real geom_area_to_m2 = config.geom_length_to_m * config.geom_length_to_m;
    Real geom_volume_to_m3 = geom_area_to_m2 * config.geom_length_to_m;
    state.volume_total_si = state.volume_total_geom * geom_volume_to_m3;
    state.projection_length_si = state.projection_length_geom * config.geom_length_to_m;

    state.effective_area_used_geom = config.effective_area_input_geom;
    if (state.effective_area_used_geom <= TinyReal)
    {
        state.effective_area_used_geom =
            state.volume_total_geom / (state.projection_length_geom + TinyReal);
    }
    state.effective_area_used_si = state.effective_area_used_geom * geom_area_to_m2;

    if (config.use_voltage_driven_source && state.effective_area_used_si > TinyReal)
    {
        state.used_voltage_driven_source = true;
        Real resistance_total =
            SMAX(static_cast<Real>(0.0),
                 config.circuit_resistance_ohm + config.external_series_resistance_ohm);
        Real inductance_total =
            SMAX(static_cast<Real>(0.0),
                 config.circuit_inductance_h + config.external_series_inductance_h);
        Real inductive_reactance =
            SMAX(static_cast<Real>(0.0), angular_frequency) * inductance_total;
        Real impedance_magnitude =
            sqrt(resistance_total * resistance_total +
                 inductive_reactance * inductive_reactance);

        state.circuit_total_resistance_ohm = resistance_total;
        state.circuit_total_inductive_reactance_ohm = inductive_reactance;
        state.circuit_impedance_magnitude_ohm = impedance_magnitude;

        if (config.voltage_rms > TinyReal && impedance_magnitude > TinyReal)
        {
            Real current_rms = config.voltage_rms / impedance_magnitude;
            Real current_peak = config.current_peak_factor * current_rms;
            Real current_phase = -atan2(inductive_reactance, resistance_total + TinyReal);
            Real source_current_peak_density =
                config.turns * current_peak /
                (state.effective_area_used_si + TinyReal);
            Real source_current_real_density = source_current_peak_density * cos(current_phase);
            Real source_current_imag_density = source_current_peak_density * sin(current_phase);

            state.frequency_source_real =
                source_current_real_density * state.source_direction_used;
            state.frequency_source_imag =
                source_current_imag_density * state.source_direction_used;
            state.harmonic_source_amplitude_used =
                source_current_peak_density * state.source_direction_used;
            state.circuit_current_rms = current_rms;
            state.circuit_current_peak = current_peak;
            state.circuit_current_phase_rad = current_phase;
        }
    }
    else if (config.use_current_driven_source && state.effective_area_used_si > TinyReal)
    {
        state.used_current_driven_source = true;
        if (config.current_rms > TinyReal)
        {
            Real current_peak = config.current_peak_factor * config.current_rms;
            Real source_current_peak_density =
                config.turns * current_peak /
                (state.effective_area_used_si + TinyReal);
            state.frequency_source_real =
                source_current_peak_density * state.source_direction_used;
            state.frequency_source_imag =
                config.imag_scale * state.frequency_source_real;
            state.harmonic_source_amplitude_used = state.frequency_source_real;
            state.circuit_current_rms = config.current_rms;
            state.circuit_current_peak = current_peak;
            state.circuit_current_phase_rad = atan2(config.imag_scale, 1.0);
        }
    }

    state.harmonic_source_magnitude_used = state.harmonic_source_amplitude_used.norm();
    state.frequency_source_real_magnitude = state.frequency_source_real.norm();
    state.frequency_source_imag_magnitude = state.frequency_source_imag.norm();
    return state;
}

inline CircularHarmonicSourceCurrentDensity::
    CircularHarmonicSourceCurrentDensity(SPHBody &sph_body,
                                        Real amplitude_magnitude,
                                        Real frequency_hz,
                                        Real &physical_time,
                                        const Vecd &coil_center,
                                        const Vecd &coil_axis,
                                        Real phase)
    : LocalDynamics(sph_body),
      amplitude_magnitude_(amplitude_magnitude),
      frequency_hz_(frequency_hz),
      phase_(phase),
      physical_time_(physical_time),
      coil_center_(coil_center),
      coil_axis_(NormalizeVectorOrDefault(coil_axis, Vecd(0.0, 0.0, 1.0))),
      positions_(particles_->getVariableDataByName<Vecd>("Position")),
      source_current_density_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensity")) {}

inline void CircularHarmonicSourceCurrentDensity::update(size_t index_i, Real dt)
{
    (void)dt;
    Real harmonic_factor = sin(2.0 * Pi * frequency_hz_ * physical_time_ + phase_);
    source_current_density_[index_i] =
        amplitude_magnitude_ * harmonic_factor *
        TangentialDirectionAroundAxis(positions_[index_i], coil_center_, coil_axis_);
}

inline CircularFrequencySourceCurrentDensity::
    CircularFrequencySourceCurrentDensity(SPHBody &sph_body,
                                         Real real_magnitude,
                                         Real imag_magnitude,
                                         const Vecd &coil_center,
                                         const Vecd &coil_axis)
    : LocalDynamics(sph_body),
      real_magnitude_(real_magnitude),
      imag_magnitude_(imag_magnitude),
      coil_center_(coil_center),
      coil_axis_(NormalizeVectorOrDefault(coil_axis, Vecd(0.0, 0.0, 1.0))),
      positions_(particles_->getVariableDataByName<Vecd>("Position")),
      source_current_density_real_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityReal")),
      source_current_density_imag_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensityImag")) {}

inline void CircularFrequencySourceCurrentDensity::update(size_t index_i, Real dt)
{
    (void)dt;
    Vecd tangent =
        TangentialDirectionAroundAxis(positions_[index_i], coil_center_, coil_axis_);
    source_current_density_real_[index_i] = real_magnitude_ * tangent;
    source_current_density_imag_[index_i] = imag_magnitude_ * tangent;
}
} // namespace extra_electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_MULTITURN_COIL_DRIVE_HPP
