/**
 * @file rh200_excitation_joule_audit.h
 * @brief Task A/C: excitation-to-Joule and heating-rate audit CSV helpers for RH200.
 */
#ifndef RH200_EXCITATION_JOULE_AUDIT_H
#define RH200_EXCITATION_JOULE_AUDIT_H

#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_observables.h"
#include "rh200_fake_joule_heat.h"
#include "rh200_joule_heat_grid.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace SPH
{
namespace rh200
{

inline constexpr const char *kExcitationToJouleAuditCsv = "./output/RH200_EXCITATION_TO_JOULE_AUDIT.csv";
inline constexpr const char *kHeatingRateAuditCsv = "./output/RH200_HEATING_RATE_AUDIT.csv";
inline constexpr Real kDefaultTargetTemperatureForEstimateK = Real(1500);

struct Rh200ScalarFieldStats
{
    Real min_val = Real(0);
    Real max_val = Real(0);
    Real mean = Real(0);
    Real vol_weighted_rms = Real(0);
};

struct Rh200VecdFieldStats
{
    Real vol_weighted_rms = Real(0);
    Real max_magnitude = Real(0);
};

struct Rh200GlassEmFieldAuditSnapshot
{
    Rh200VecdFieldStats a_coil_real_on_glass;
    Rh200VecdFieldStats a_coil_imag_on_glass;
    Rh200ScalarFieldStats phi_real;
    Rh200ScalarFieldStats phi_imag;
    Rh200VecdFieldStats e_real;
    Rh200VecdFieldStats e_imag;
    Rh200VecdFieldStats j_real;
    Rh200VecdFieldStats j_imag;
    Rh200ScalarFieldStats joule_heat;
    Real sigma_min = Real(0);
    Real sigma_max = Real(0);
    Real sigma_mean = Real(0);
    Real p_joule_particle = Real(0);
};

struct Rh200GlassThermalMassAudit
{
    Real glass_volume = Real(0);
    Real glass_mass = Real(0);
    Real cp = Real(0);
    Real mcp = Real(0);
};

struct Rh200ExcitationToJouleAuditRecord
{
    std::string case_name = "rh200_glass_em_stirring";
    std::string joule_mode = "off";
    std::string em_mode = "edge-flux-complex";
    Real frequency_hz = Real(0);
    Real omega = Real(0);
    Real coil_radius = Real(0);
    Real coil_z_min = Real(0);
    Real coil_z_max = Real(0);
    size_t coil_turns = 0;
    Real coil_center_x = Real(0);
    Real coil_center_y = Real(0);
    Real coil_center_z = Real(0);
    Real base_current_per_loop = Real(0);
    Real equivalent_current_scale = Real(1);
    Real em_power_scale_factor = Real(1);
    Real grid_rescale_factor = Real(1);
    Real target_power = Real(0);
    Rh200GlassEmFieldAuditSnapshot raw_fields;
    Rh200GlassEmFieldAuditSnapshot scaled_fields;
    Real p_joule_raw_grid_sample_initial = Real(0);
    Real p_joule_scaled_grid_sample_initial = Real(0);
    Rh200GlassThermalMassAudit thermal_mass;
    Real stl_glass_volume = Real(0);
};

struct Rh200HeatingRateAuditRecord
{
    Real time = Real(0);
    Real glass_volume = Real(0);
    Real glass_mass = Real(0);
    Real cp = Real(0);
    Real mcp = Real(0);
    Real target_power = Real(0);
    Real p_grid_sample_current = Real(0);
    Real integrated_joule_energy = Real(0);
    Real thermal_energy = Real(0);
    Real t_mean = Real(0);
    Real t_min = Real(0);
    Real t_max = Real(0);
    Real t_mean_initial = Real(0);
    Real dT_mean_measured = Real(0);
    Real dT_mean_expected = Real(0);
    Real heating_rate_expected = Real(0);
    Real target_temperature_for_estimate_K = kDefaultTargetTemperatureForEstimateK;
    Real time_to_target_temperature_no_loss_s = Real(0);
    Real out_of_grid_particle_count = Real(0);
    Real out_of_grid_particle_fraction = Real(0);
};

inline Real hostVolWeightedVecdNormSquaredLocal(BaseParticles &particles, const std::string &field_name, size_t n)
{
    electromagnetics::ophelie::syncVariableToHost<Vecd>(particles, field_name);
    electromagnetics::ophelie::syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *values = particles.getVariableDataByName<Vecd>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        sum += vol[i] * values[i].squaredNorm();
    }
    return sum;
}

inline Real hostVecdVolWeightedNormLocal(BaseParticles &particles, const std::string &field_name, size_t n)
{
    return std::sqrt(std::max(hostVolWeightedVecdNormSquaredLocal(particles, field_name, n), Real(0)));
}

inline Rh200VecdFieldStats measureRh200VecdFieldStats(BaseParticles &particles, const std::string &field_name, size_t n)
{
    Rh200VecdFieldStats stats;
    stats.vol_weighted_rms = hostVecdVolWeightedNormLocal(particles, field_name, n);
    stats.max_magnitude = electromagnetics::ophelie::hostVecdFieldMax(particles, field_name, n);
    return stats;
}

inline Rh200ScalarFieldStats measureRh200ScalarFieldStats(BaseParticles &particles, const std::string &field_name,
                                                          size_t n)
{
    Rh200ScalarFieldStats stats;
    electromagnetics::ophelie::syncVariableToHost<Real>(particles, field_name);
    electromagnetics::ophelie::syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *values = particles.getVariableDataByName<Real>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    if (n == 0)
    {
        return stats;
    }
    Real vol_sum = Real(0);
    Real weighted_sum = Real(0);
    stats.min_val = values[0];
    stats.max_val = values[0];
    for (size_t i = 0; i < n; ++i)
    {
        stats.min_val = std::min(stats.min_val, values[i]);
        stats.max_val = std::max(stats.max_val, values[i]);
        vol_sum += vol[i];
        weighted_sum += vol[i] * values[i];
    }
    stats.mean = weighted_sum / (vol_sum + TinyReal);
    Real sum_sq = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        sum_sq += vol[i] * values[i] * values[i];
    }
    stats.vol_weighted_rms = std::sqrt(std::max(sum_sq, Real(0)));
    return stats;
}

inline void measureRh200SigmaStats(BaseParticles &particles, const std::string &sigma_field, size_t n, Real &sigma_min,
                                   Real &sigma_max, Real &sigma_mean)
{
    const Rh200ScalarFieldStats stats = measureRh200ScalarFieldStats(particles, sigma_field, n);
    sigma_min = stats.min_val;
    sigma_max = stats.max_val;
    sigma_mean = stats.mean;
}

inline Rh200GlassEmFieldAuditSnapshot collectRh200GlassEmFieldAuditSnapshot(
    BaseParticles &particles, const electromagnetics::ophelie::OphelieGlassFieldNames &glass_names, size_t n)
{
    Rh200GlassEmFieldAuditSnapshot snap;
    snap.a_coil_real_on_glass = measureRh200VecdFieldStats(particles, glass_names.a_coil_real, n);
    snap.a_coil_imag_on_glass = measureRh200VecdFieldStats(particles, glass_names.a_coil_imag, n);
    snap.phi_real = measureRh200ScalarFieldStats(particles, glass_names.phi_real, n);
    snap.phi_imag = measureRh200ScalarFieldStats(particles, glass_names.phi_imag, n);
    snap.e_real = measureRh200VecdFieldStats(particles, glass_names.e_real, n);
    snap.e_imag = measureRh200VecdFieldStats(particles, glass_names.e_imag, n);
    snap.j_real = measureRh200VecdFieldStats(particles, glass_names.j_real, n);
    snap.j_imag = measureRh200VecdFieldStats(particles, glass_names.j_imag, n);
    snap.joule_heat = measureRh200ScalarFieldStats(particles, glass_names.joule_heat, n);
    measureRh200SigmaStats(particles, glass_names.sigma, n, snap.sigma_min, snap.sigma_max, snap.sigma_mean);
    snap.p_joule_particle = electromagnetics::ophelie::hostVolWeightedSum(particles, glass_names.joule_heat, n);
    return snap;
}

inline Rh200GlassThermalMassAudit measureRh200GlassParticleThermalMass(BaseParticles &particles, Real rho, Real cp)
{
    Rh200GlassThermalMassAudit audit;
    audit.cp = cp;
    electromagnetics::ophelie::syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        audit.glass_volume += vol[i];
        audit.glass_mass += rho * vol[i];
    }
    audit.mcp = audit.glass_mass * audit.cp;
    return audit;
}

inline Real hostSampleGridPowerFromDeposition(const Rh200JouleHeatGridSpec &spec,
                                              const Rh200EmParticleDepositionHost &dep, bool apply_grid_rescale,
                                              Real grid_rescale_factor, Real &grid_rescale_out)
{
    Rh200ScalarEulerianGrid grid;
    grid.spec_ = spec;
    grid.resetAccumulators();
    grid.depositScalarCloudInCell(dep.position, dep.joule_heat, dep.volumetric_measure);
    grid.finalizeFromAccumulators();
    const Real p_before = hostSamplePowerFromScalarGrid(grid, dep.position, dep.volumetric_measure);
    grid_rescale_out = Real(1);
    if (apply_grid_rescale && p_before > TinyReal)
    {
        grid_rescale_out = grid_rescale_factor;
        grid.scaleField(grid_rescale_out);
    }
    return hostSamplePowerFromScalarGrid(grid, dep.position, dep.volumetric_measure);
}

inline Real hostSampleRawGridPowerFromDeposition(const Rh200JouleHeatGridSpec &spec,
                                                 const Rh200EmParticleDepositionHost &dep)
{
    Real unused = Real(1);
    return hostSampleGridPowerFromDeposition(spec, dep, false, Real(1), unused);
}

inline void writeRh200ExcitationToJouleAuditCsv(const std::string &path, const Rh200ExcitationToJouleAuditRecord &rec)
{
    const auto &raw = rec.raw_fields;
    const auto &scaled = rec.scaled_fields;
    const Real expected_heating_rate_raw =
        rec.thermal_mass.mcp > TinyReal ? raw.p_joule_particle / rec.thermal_mass.mcp : Real(0);
    const Real expected_heating_rate_scaled =
        rec.thermal_mass.mcp > TinyReal ? scaled.p_joule_particle / rec.thermal_mass.mcp : Real(0);

    std::ofstream out(path);
    out << std::setprecision(12);
    out << "case_name,joule_mode,em_mode,frequency_hz,omega,coil_radius,coil_z_min,coil_z_max,coil_turns,"
        << "coil_center_x,coil_center_y,coil_center_z,base_current_per_loop,equivalent_current_scale,"
        << "em_power_scale_factor,grid_rescale_factor,target_power_W,"
        << "sigma_min,sigma_max,sigma_mean,"
        << "ACoilReal_on_glass_rms_raw,ACoilReal_on_glass_max_raw,ACoilImag_on_glass_rms_raw,ACoilImag_on_glass_max_raw,"
        << "PhiReal_rms_raw,PhiImag_rms_raw,EReal_rms_raw,EReal_max_raw,EImag_rms_raw,EImag_max_raw,"
        << "JReal_rms_raw,JReal_max_raw,JImag_rms_raw,JImag_max_raw,"
        << "JouleHeat_raw_min,JouleHeat_raw_max,JouleHeat_raw_mean,P_joule_raw_particle,P_joule_raw_grid_sample_initial,"
        << "ACoilReal_on_glass_rms_scaled,ACoilReal_on_glass_max_scaled,ACoilImag_on_glass_rms_scaled,"
        << "ACoilImag_on_glass_max_scaled,PhiReal_rms_scaled,PhiImag_rms_scaled,EReal_rms_scaled,EReal_max_scaled,"
        << "EImag_rms_scaled,EImag_max_scaled,JReal_rms_scaled,JReal_max_scaled,JImag_rms_scaled,JImag_max_scaled,"
        << "JouleHeat_scaled_min,JouleHeat_scaled_max,JouleHeat_scaled_mean,P_joule_scaled_particle,"
        << "P_joule_scaled_grid_sample_initial,glass_volume_particle,glass_mass_particle,cp_J_per_kgK,Mcp_J_per_K,"
        << "stl_glass_volume_m3,expected_heating_rate_raw_K_per_s,expected_heating_rate_scaled_K_per_s\n";
    out << rec.case_name << "," << rec.joule_mode << "," << rec.em_mode << "," << rec.frequency_hz << "," << rec.omega
        << "," << rec.coil_radius << "," << rec.coil_z_min << "," << rec.coil_z_max << "," << rec.coil_turns << ","
        << rec.coil_center_x << "," << rec.coil_center_y << "," << rec.coil_center_z << "," << rec.base_current_per_loop
        << "," << rec.equivalent_current_scale << "," << rec.em_power_scale_factor << "," << rec.grid_rescale_factor
        << "," << rec.target_power << "," << raw.sigma_min << "," << raw.sigma_max << "," << raw.sigma_mean << ","
        << raw.a_coil_real_on_glass.vol_weighted_rms << "," << raw.a_coil_real_on_glass.max_magnitude << ","
        << raw.a_coil_imag_on_glass.vol_weighted_rms << "," << raw.a_coil_imag_on_glass.max_magnitude << ","
        << raw.phi_real.vol_weighted_rms << "," << raw.phi_imag.vol_weighted_rms << "," << raw.e_real.vol_weighted_rms
        << "," << raw.e_real.max_magnitude << "," << raw.e_imag.vol_weighted_rms << "," << raw.e_imag.max_magnitude
        << "," << raw.j_real.vol_weighted_rms << "," << raw.j_real.max_magnitude << "," << raw.j_imag.vol_weighted_rms
        << "," << raw.j_imag.max_magnitude << "," << raw.joule_heat.min_val << "," << raw.joule_heat.max_val << ","
        << raw.joule_heat.mean << "," << raw.p_joule_particle << "," << rec.p_joule_raw_grid_sample_initial << ","
        << scaled.a_coil_real_on_glass.vol_weighted_rms << "," << scaled.a_coil_real_on_glass.max_magnitude << ","
        << scaled.a_coil_imag_on_glass.vol_weighted_rms << "," << scaled.a_coil_imag_on_glass.max_magnitude << ","
        << scaled.phi_real.vol_weighted_rms << "," << scaled.phi_imag.vol_weighted_rms << ","
        << scaled.e_real.vol_weighted_rms << "," << scaled.e_real.max_magnitude << ","
        << scaled.e_imag.vol_weighted_rms << "," << scaled.e_imag.max_magnitude << ","
        << scaled.j_real.vol_weighted_rms << "," << scaled.j_real.max_magnitude << ","
        << scaled.j_imag.vol_weighted_rms << "," << scaled.j_imag.max_magnitude << ","
        << scaled.joule_heat.min_val << "," << scaled.joule_heat.max_val << "," << scaled.joule_heat.mean << ","
        << scaled.p_joule_particle << "," << rec.p_joule_scaled_grid_sample_initial << ","
        << rec.thermal_mass.glass_volume << "," << rec.thermal_mass.glass_mass << "," << rec.thermal_mass.cp << ","
        << rec.thermal_mass.mcp << "," << rec.stl_glass_volume << "," << expected_heating_rate_raw << ","
        << expected_heating_rate_scaled << "\n";
    std::cout << "[rh200] wrote excitation-to-Joule audit: " << path << std::endl;
}

inline void writeRh200HeatingRateAuditCsvHeader(const std::string &path)
{
    std::ofstream out(path, std::ios::out);
    out << "time,glass_volume,glass_mass,cp,Mcp,target_power,P_grid_sample_current,integrated_joule_energy,"
        << "thermal_energy,T_mean,T_min,T_max,T_mean_initial,dT_mean_measured,dT_mean_expected,heating_rate_expected,"
        << "target_temperature_for_estimate_K,time_to_target_temperature_no_loss_s,out_of_grid_particle_count,"
        << "out_of_grid_particle_fraction\n";
}

inline Real hostGlassTemperatureMin(BaseParticles &particles, const std::string &temperature_name)
{
    electromagnetics::ophelie::syncVariableToHost<Real>(particles, temperature_name);
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_name);
    const size_t n = particles.TotalRealParticles();
    if (n == 0)
    {
        return Real(0);
    }
    Real t_min = temperature[0];
    for (size_t i = 1; i < n; ++i)
    {
        t_min = std::min(t_min, temperature[i]);
    }
    return t_min;
}

inline Rh200HeatingRateAuditRecord makeRh200HeatingRateAuditRecord(
    Real time, const Rh200GlassThermalMassAudit &thermal_mass, Real target_power, Real p_grid_sample_current,
    Real integrated_joule_energy, Real thermal_energy, Real t_mean, Real t_min, Real t_max, Real t_mean_initial,
    Real out_of_grid_count, Real out_of_grid_fraction,
    Real target_temperature_for_estimate_K = kDefaultTargetTemperatureForEstimateK)
{
    Rh200HeatingRateAuditRecord rec;
    rec.time = time;
    rec.glass_volume = thermal_mass.glass_volume;
    rec.glass_mass = thermal_mass.glass_mass;
    rec.cp = thermal_mass.cp;
    rec.mcp = thermal_mass.mcp;
    rec.target_power = target_power;
    rec.p_grid_sample_current = p_grid_sample_current;
    rec.integrated_joule_energy = integrated_joule_energy;
    rec.thermal_energy = thermal_energy;
    rec.t_mean = t_mean;
    rec.t_min = t_min;
    rec.t_max = t_max;
    rec.t_mean_initial = t_mean_initial;
    rec.dT_mean_measured = t_mean - t_mean_initial;
    rec.dT_mean_expected = thermal_mass.mcp > TinyReal ? integrated_joule_energy / thermal_mass.mcp : Real(0);
    rec.heating_rate_expected = thermal_mass.mcp > TinyReal ? p_grid_sample_current / thermal_mass.mcp : Real(0);
    rec.target_temperature_for_estimate_K = target_temperature_for_estimate_K;
    if (p_grid_sample_current > TinyReal && target_temperature_for_estimate_K > t_mean)
    {
        rec.time_to_target_temperature_no_loss_s =
            thermal_mass.mcp * (target_temperature_for_estimate_K - t_mean) / p_grid_sample_current;
    }
    else
    {
        rec.time_to_target_temperature_no_loss_s = Real(0);
    }
    rec.out_of_grid_particle_count = out_of_grid_count;
    rec.out_of_grid_particle_fraction = out_of_grid_fraction;
    return rec;
}

inline void appendRh200HeatingRateAuditCsv(const std::string &path, const Rh200HeatingRateAuditRecord &rec)
{
    std::ofstream out(path, std::ios::app);
    out << std::setprecision(12);
    out << rec.time << "," << rec.glass_volume << "," << rec.glass_mass << "," << rec.cp << "," << rec.mcp << ","
        << rec.target_power << "," << rec.p_grid_sample_current << "," << rec.integrated_joule_energy << ","
        << rec.thermal_energy << "," << rec.t_mean << "," << rec.t_min << "," << rec.t_max << "," << rec.t_mean_initial
        << "," << rec.dT_mean_measured << "," << rec.dT_mean_expected << "," << rec.heating_rate_expected << ","
        << rec.target_temperature_for_estimate_K << "," << rec.time_to_target_temperature_no_loss_s << ","
        << rec.out_of_grid_particle_count << "," << rec.out_of_grid_particle_fraction << "\n";
}

} // namespace rh200
} // namespace SPH

#endif // RH200_EXCITATION_JOULE_AUDIT_H
