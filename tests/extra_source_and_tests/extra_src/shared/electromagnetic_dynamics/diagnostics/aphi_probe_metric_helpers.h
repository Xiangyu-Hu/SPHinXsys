#ifndef APHI_PROBE_METRIC_HELPERS_H
#define APHI_PROBE_METRIC_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_team7_canonical_case_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiProbeLineSample
{
    Vecd position = Vecd::Zero();
    size_t nearest_particle_index = 0;
    Real neighbor_distance = 0.0;
    Vecd b_real = Vecd::Zero();
    Vecd b_imag = Vecd::Zero();
    Vecd j_real = Vecd::Zero();
    Vecd j_imag = Vecd::Zero();
    Real b_abs = 0.0;
    Real j_abs = 0.0;
    Real joule = 0.0;
    Real sigma = 0.0;
    Real region_id = 0.0;
};

struct AphiProbeRegionMetricsRow
{
    std::string region_name;
    AphiBenchmarkMaterialRegion region = AphiBenchmarkMaterialRegion::AllCore;
    size_t particle_count = 0;
    Real joule_integral = 0.0;
    Real b_max = 0.0;
    Real b_avg = 0.0;
    Real j_max = 0.0;
    Real j_avg = 0.0;
    Real joule_max = 0.0;
    Real joule_avg = 0.0;
    Real e_max = 0.0;
    Real sigma_avg = 0.0;
};

struct AphiProbeObservedFieldNames
{
    std::string b_real = "BReal";
    std::string b_imag = "BImag";
    std::string b_magnitude = "BMagnitude";
    std::string e_magnitude = "EMagnitude";
    std::string j_magnitude = "JMagnitude";
    std::string material_region_id = "MaterialRegionId";
};

inline bool writeGmresSummaryCsv(const std::string &path, bool converged, UnsignedInt num_iterations, Real final_residual,
                                 Real omega, Real dp, Real body_length, Real body_height, Real body_width,
                                 size_t particles, Real conductor_joule_integral, Real air_joule_integral,
                                 Real total_joule_integral)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "converged,num_iterations,final_residual,omega,dp,body_length,body_height,body_width,particles,"
              "conductor_Joule_integral,air_Joule_integral,total_Joule_integral\n";
    output << (converged ? 1 : 0) << "," << num_iterations << "," << final_residual << "," << omega << "," << dp
           << "," << body_length << "," << body_height << "," << body_width << "," << particles << ","
           << conductor_joule_integral << "," << air_joule_integral << "," << total_joule_integral << "\n";
    return true;
}

struct AphiProbeCsvWriteResult
{
    std::string output_dir = "output";
    size_t b_line_samples = 0;
    size_t j_line_samples = 0;
    size_t point_probe_count = 0;
    size_t region_rows = 0;
    bool gmres_summary_written = false;
    bool all_probe_values_finite = false;
};

inline Real vec3Abs(const Vecd &value) { return value.norm(); }

inline size_t nearestRealParticleIndex(const Vecd *positions, size_t total_real_particles, const Vecd &target)
{
    size_t nearest_index = 0;
    Real nearest_distance_squared = std::numeric_limits<Real>::max();
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        const Real distance_squared = (positions[i] - target).squaredNorm();
        if (distance_squared < nearest_distance_squared)
        {
            nearest_distance_squared = distance_squared;
            nearest_index = i;
        }
    }
    return nearest_index;
}

inline StdVec<AphiProbeLineSample> sampleLineProbesNearestParticle(
    BaseParticles &particles, const Vecd *positions, size_t total_real_particles, const StdVec<Vecd> &sample_positions,
    const AphiProbeObservedFieldNames &obs_fields, const AphiJouleHeatingFieldNames &joule_fields,
    const AphiMaterialNames &material_names)
{
    syncVariableToHost<Vecd>(particles, obs_fields.b_real);
    syncVariableToHost<Vecd>(particles, obs_fields.b_imag);
    syncVariableToHost<Vecd>(particles, joule_fields.current_density_real);
    syncVariableToHost<Vecd>(particles, joule_fields.current_density_imag);
    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    syncVariableToHost<Real>(particles, material_names.sigma);
    syncVariableToHost<Real>(particles, obs_fields.material_region_id);

    const Vecd *b_real = particles.getVariableDataByName<Vecd>(obs_fields.b_real);
    const Vecd *b_imag = particles.getVariableDataByName<Vecd>(obs_fields.b_imag);
    const Vecd *j_real = particles.getVariableDataByName<Vecd>(joule_fields.current_density_real);
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(joule_fields.current_density_imag);
    const Real *joule = particles.getVariableDataByName<Real>(joule_fields.joule_heat_source);
    const Real *sigma = particles.getVariableDataByName<Real>(material_names.sigma);
    const Real *region_id = particles.getVariableDataByName<Real>(obs_fields.material_region_id);

    StdVec<AphiProbeLineSample> samples;
    samples.reserve(sample_positions.size());
    for (const Vecd &sample_position : sample_positions)
    {
        const size_t index = nearestRealParticleIndex(positions, total_real_particles, sample_position);
        AphiProbeLineSample sample;
        sample.position = sample_position;
        sample.nearest_particle_index = index;
        sample.neighbor_distance = (positions[index] - sample_position).norm();
        sample.b_real = b_real[index];
        sample.b_imag = b_imag[index];
        sample.j_real = j_real[index];
        sample.j_imag = j_imag[index];
        sample.b_abs = vec3Abs(sample.b_real);
        sample.j_abs = vec3Abs(sample.j_real);
        sample.joule = joule[index];
        sample.sigma = sigma[index];
        sample.region_id = region_id[index];
        samples.push_back(sample);
    }
    return samples;
}

inline StdVec<Vecd> buildTeam7CenterlineSamplePositions(Real body_length, Real body_height, Real body_width,
                                                        UnsignedInt sample_count)
{
    const Vecd center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    StdVec<Vecd> positions;
    positions.reserve(static_cast<size_t>(sample_count));
    for (UnsignedInt sample_index = 0; sample_index < sample_count; ++sample_index)
    {
        const Real alpha =
            static_cast<Real>(sample_index) / (static_cast<Real>(sample_count - 1) + TinyReal);
        positions.push_back(Vecd(alpha * body_length, center[1], center[2]));
    }
    return positions;
}

inline bool probeSamplesFinite(const StdVec<AphiProbeLineSample> &samples)
{
    for (const AphiProbeLineSample &sample : samples)
    {
        if (!sample.b_real.allFinite() || !sample.b_imag.allFinite() || !sample.j_real.allFinite() ||
            !sample.j_imag.allFinite() || !std::isfinite(sample.b_abs) || !std::isfinite(sample.j_abs) ||
            !std::isfinite(sample.joule) || !std::isfinite(sample.sigma) || !std::isfinite(sample.region_id))
        {
            return false;
        }
    }
    return !samples.empty();
}

inline AphiProbeRegionMetricsRow hostRegionProbeMetricsRow(
    BaseParticles &particles, const Vecd *positions, size_t total_real_particles,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, AphiBenchmarkMaterialRegion region,
    const std::string &region_name, const AphiProbeObservedFieldNames &obs_fields,
    const AphiJouleHeatingFieldNames &joule_fields, const AphiMaterialNames &material_names)
{
    const auto in_region = [&](const Vecd &position) {
        return team7ParticleInRegion(position, layout, region);
    };

    AphiProbeRegionMetricsRow row;
    row.region_name = region_name;
    row.region = region;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_region(positions[i]))
        {
            ++row.particle_count;
        }
    }
    row.joule_integral = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles, in_region,
                                                                 joule_fields.joule_heat_source);
    row.b_max = hostParticleRegionScalarMax(particles, obs_fields.b_magnitude, positions, total_real_particles, in_region);
    row.j_max = hostParticleRegionScalarMax(particles, obs_fields.j_magnitude, positions, total_real_particles, in_region);
    row.joule_max = hostParticleRegionScalarMax(particles, joule_fields.joule_heat_source, positions, total_real_particles,
                                                in_region);
    row.e_max = hostParticleRegionScalarMax(particles, obs_fields.e_magnitude, positions, total_real_particles, in_region);

    syncVariableToHost<Real>(particles, obs_fields.b_magnitude);
    syncVariableToHost<Real>(particles, obs_fields.j_magnitude);
    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    syncVariableToHost<Real>(particles, material_names.sigma);
    const Real *b_mag = particles.getVariableDataByName<Real>(obs_fields.b_magnitude);
    const Real *j_mag = particles.getVariableDataByName<Real>(obs_fields.j_magnitude);
    const Real *joule = particles.getVariableDataByName<Real>(joule_fields.joule_heat_source);
    const Real *sigma = particles.getVariableDataByName<Real>(material_names.sigma);

    Real b_sum = 0.0;
    Real j_sum = 0.0;
    Real joule_sum = 0.0;
    Real sigma_sum = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_region(positions[i]))
        {
            continue;
        }
        b_sum += b_mag[i];
        j_sum += j_mag[i];
        joule_sum += joule[i];
        sigma_sum += sigma[i];
    }
    if (row.particle_count > 0)
    {
        const Real count = static_cast<Real>(row.particle_count);
        row.b_avg = b_sum / count;
        row.j_avg = j_sum / count;
        row.joule_avg = joule_sum / count;
        row.sigma_avg = sigma_sum / count;
    }
    return row;
}

inline bool writeProbeLineCsv(const std::string &path, const StdVec<AphiProbeLineSample> &samples, bool include_j)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    if (include_j)
    {
        output << "x,y,z,nearest_index,nearest_distance,"
               << "Bx_real,By_real,Bz_real,Bx_imag,By_imag,Bz_imag,B_abs,"
               << "Jx_real,Jy_real,Jz_real,Jx_imag,Jy_imag,Jz_imag,J_abs,Joule,Sigma,RegionId\n";
        for (const AphiProbeLineSample &sample : samples)
        {
            output << sample.position[0] << "," << sample.position[1] << "," << sample.position[2] << ","
                   << sample.nearest_particle_index << "," << sample.neighbor_distance << "," << sample.b_real[0] << ","
                   << sample.b_real[1] << "," << sample.b_real[2] << "," << sample.b_imag[0] << "," << sample.b_imag[1]
                   << "," << sample.b_imag[2] << "," << sample.b_abs << "," << sample.j_real[0] << ","
                   << sample.j_real[1] << "," << sample.j_real[2] << "," << sample.j_imag[0] << "," << sample.j_imag[1]
                   << "," << sample.j_imag[2] << "," << sample.j_abs << "," << sample.joule << "," << sample.sigma
                   << "," << sample.region_id << "\n";
        }
    }
    else
    {
        output << "x,y,z,nearest_index,nearest_distance,"
               << "Bx_real,By_real,Bz_real,Bx_imag,By_imag,Bz_imag,B_abs,Joule,Sigma,RegionId\n";
        for (const AphiProbeLineSample &sample : samples)
        {
            output << sample.position[0] << "," << sample.position[1] << "," << sample.position[2] << ","
                   << sample.nearest_particle_index << "," << sample.neighbor_distance << "," << sample.b_real[0] << ","
                   << sample.b_real[1] << "," << sample.b_real[2] << "," << sample.b_imag[0] << "," << sample.b_imag[1]
                   << "," << sample.b_imag[2] << "," << sample.b_abs << "," << sample.joule << "," << sample.sigma
                   << "," << sample.region_id << "\n";
        }
    }
    return true;
}

inline bool writeRegionMetricsCsv(const std::string &path, const StdVec<AphiProbeRegionMetricsRow> &rows)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "region,particle_count,joule_integral,b_max,b_avg,j_max,j_avg,joule_max,joule_avg,e_max,sigma_avg\n";
    for (const AphiProbeRegionMetricsRow &row : rows)
    {
        output << row.region_name << "," << row.particle_count << "," << row.joule_integral << "," << row.b_max << ","
               << row.b_avg << "," << row.j_max << "," << row.j_avg << "," << row.joule_max << "," << row.joule_avg
               << "," << row.e_max << "," << row.sigma_avg << "\n";
    }
    return true;
}

inline AphiProbeCsvWriteResult writeSourceDrivenProbeCsvArtifacts(
    BaseParticles &particles, const benchmark::AphiTeam7LikeUnitBoxLayout &layout, Real body_length, Real body_height,
    Real body_width, const AphiProbeObservedFieldNames &obs_fields, const AphiJouleHeatingFieldNames &joule_fields,
    const AphiMaterialNames &material_names, bool converged, UnsignedInt num_iterations, Real final_residual, Real omega,
    Real dp, size_t particle_count, Real conductor_joule_integral, Real air_joule_integral, Real total_joule_integral,
    const std::string &output_dir = "output", UnsignedInt centerline_sample_count = 32)
{
    AphiProbeCsvWriteResult result;
    result.output_dir = output_dir;

    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    const StdVec<Vecd> centerline_positions =
        buildTeam7CenterlineSamplePositions(body_length, body_height, body_width, centerline_sample_count);
    const StdVec<AphiProbeLineSample> line_samples = sampleLineProbesNearestParticle(
        particles, positions, total_real_particles, centerline_positions, obs_fields, joule_fields, material_names);
    result.b_line_samples = line_samples.size();
    result.j_line_samples = line_samples.size();
    result.all_probe_values_finite = probeSamplesFinite(line_samples);

    StdVec<AphiProbeRegionMetricsRow> region_rows;
    region_rows.push_back(hostRegionProbeMetricsRow(particles, positions, total_real_particles, layout,
                                                    AphiBenchmarkMaterialRegion::Air, "air", obs_fields, joule_fields,
                                                    material_names));
    region_rows.push_back(hostRegionProbeMetricsRow(particles, positions, total_real_particles, layout,
                                                    AphiBenchmarkMaterialRegion::Coil, "coil", obs_fields, joule_fields,
                                                    material_names));
    region_rows.push_back(hostRegionProbeMetricsRow(particles, positions, total_real_particles, layout,
                                                    AphiBenchmarkMaterialRegion::Conductor, "conductor", obs_fields,
                                                    joule_fields, material_names));
    region_rows.push_back(hostRegionProbeMetricsRow(particles, positions, total_real_particles, layout,
                                                    AphiBenchmarkMaterialRegion::AllCore, "all", obs_fields, joule_fields,
                                                    material_names));
    result.region_rows = region_rows.size();

    const bool b_written = writeProbeLineCsv(output_dir + "/aphi_probe_B_line.csv", line_samples, false);
    const bool j_written = writeProbeLineCsv(output_dir + "/aphi_probe_J_line.csv", line_samples, true);
    const bool region_written = writeRegionMetricsCsv(output_dir + "/aphi_region_metrics.csv", region_rows);
    result.gmres_summary_written =
        writeGmresSummaryCsv(output_dir + "/aphi_gmres_summary.csv", converged, num_iterations, final_residual, omega, dp,
                               body_length, body_height, body_width, particle_count, conductor_joule_integral,
                               air_joule_integral, total_joule_integral) &&
        b_written && j_written && region_written;
    result.point_probe_count = line_samples.size();
    return result;
}

inline bool probeMetricCsvPassed(const AphiProbeCsvWriteResult &result)
{
    return result.gmres_summary_written && result.b_line_samples > 0 && result.j_line_samples > 0 &&
           result.region_rows >= 3 && result.all_probe_values_finite;
}

inline void printProbeCsvWriteResult(const char *test_name, const AphiProbeCsvWriteResult &result, bool passed)
{
    std::cout << test_name << " probe_csv_dir=" << result.output_dir << " b_line_samples=" << result.b_line_samples
              << " j_line_samples=" << result.j_line_samples << " region_rows=" << result.region_rows
              << " all_finite=" << (result.all_probe_values_finite ? 1 : 0)
              << " gmres_summary_written=" << (result.gmres_summary_written ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
