#ifndef ELECTROMAGNETIC_MULTITURN_COIL_DRIVE_H
#define ELECTROMAGNETIC_MULTITURN_COIL_DRIVE_H

#include "sphinxsys.h"

namespace SPH
{
namespace extra_electromagnetics
{
struct MultiturnCoilDriveConfig
{
    bool use_current_driven_source = false;
    bool use_voltage_driven_source = false;
    bool use_circular_source = false;
    Real current_rms = 0.0;
    Real voltage_rms = 0.0;
    Real turns = 1.0;
    Real current_peak_factor = 1.0;
    Real imag_scale = 0.0;
    Real circuit_resistance_ohm = 0.0;
    Real circuit_inductance_h = 0.0;
    Real external_series_resistance_ohm = 0.0;
    Real external_series_inductance_h = 0.0;
    Real effective_area_input_geom = -1.0;
    Real geom_length_to_m = 1.0e-3;
    Vecd source_direction = ZeroData<Vecd>::value;
    Vecd source_axis = ZeroData<Vecd>::value;
    Vecd fallback_harmonic_amplitude = ZeroData<Vecd>::value;
};

struct MultiturnCoilDriveState
{
    bool used_current_driven_source = false;
    bool used_voltage_driven_source = false;
    Vecd source_direction_used = ZeroData<Vecd>::value;
    Vecd source_axis_used = ZeroData<Vecd>::value;
    Vecd harmonic_source_amplitude_used = ZeroData<Vecd>::value;
    Vecd frequency_source_real = ZeroData<Vecd>::value;
    Vecd frequency_source_imag = ZeroData<Vecd>::value;
    Real harmonic_source_magnitude_used = 0.0;
    Real frequency_source_real_magnitude = 0.0;
    Real frequency_source_imag_magnitude = 0.0;
    Real effective_area_used_geom = 0.0;
    Real effective_area_used_si = 0.0;
    Real volume_total_geom = 0.0;
    Real volume_total_si = 0.0;
    Real projection_length_geom = 0.0;
    Real projection_length_si = 0.0;
    Real mean_radius_geom = 0.0;
    Real circuit_total_resistance_ohm = 0.0;
    Real circuit_total_inductive_reactance_ohm = 0.0;
    Real circuit_impedance_magnitude_ohm = 0.0;
    Real circuit_current_rms = 0.0;
    Real circuit_current_peak = 0.0;
    Real circuit_current_phase_rad = 0.0;
};

Vecd NormalizeVectorOrDefault(const Vecd &candidate, const Vecd &fallback);
Vecd TangentialDirectionAroundAxis(const Vecd &position, const Vecd &center, const Vecd &axis);

MultiturnCoilDriveState BuildMultiturnCoilDriveState(BaseParticles &coil_particles,
                                                     const BoundingBoxd &coil_bounds,
                                                     const Vecd &coil_center,
                                                     Real angular_frequency,
                                                     const MultiturnCoilDriveConfig &config);

class CircularHarmonicSourceCurrentDensity : public LocalDynamics
{
  public:
    explicit CircularHarmonicSourceCurrentDensity(SPHBody &sph_body,
                                                  Real amplitude_magnitude,
                                                  Real frequency_hz,
                                                  Real &physical_time,
                                                  const Vecd &coil_center,
                                                  const Vecd &coil_axis,
                                                  Real phase = 0.0);
    virtual ~CircularHarmonicSourceCurrentDensity() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real amplitude_magnitude_;
    Real frequency_hz_;
    Real phase_;
    Real &physical_time_;
    Vecd coil_center_;
    Vecd coil_axis_;
    Vecd *positions_;
    Vecd *source_current_density_;
};

class CircularFrequencySourceCurrentDensity : public LocalDynamics
{
  public:
    explicit CircularFrequencySourceCurrentDensity(SPHBody &sph_body,
                                                   Real real_magnitude,
                                                   Real imag_magnitude,
                                                   const Vecd &coil_center,
                                                   const Vecd &coil_axis);
    virtual ~CircularFrequencySourceCurrentDensity() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real real_magnitude_;
    Real imag_magnitude_;
    Vecd coil_center_;
    Vecd coil_axis_;
    Vecd *positions_;
    Vecd *source_current_density_real_;
    Vecd *source_current_density_imag_;
};
} // namespace extra_electromagnetics
} // namespace SPH

#include "electromagnetic_multiturn_coil_drive.hpp"

#endif // ELECTROMAGNETIC_MULTITURN_COIL_DRIVE_H
