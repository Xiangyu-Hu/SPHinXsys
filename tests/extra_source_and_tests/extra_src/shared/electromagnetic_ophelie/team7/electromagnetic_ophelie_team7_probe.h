#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_PROBE_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_PROBE_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_team7_native_geometry.h"

#include <cmath>
#include <filesystem>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

namespace fs = std::filesystem;

enum class Team7BzProbeLine
{
    A1B1,
    A2B2
};

struct Team7BzProbeLineSpec
{
    Real y_mm = 72.0;
    Real z_mm = 34.0;
    const char *csv_basename = "TEAM7_Bz_A1_B1_reference_mT.csv";
};

inline Team7BzProbeLineSpec team7BzProbeLineSpec(Team7BzProbeLine line)
{
    Team7BzProbeLineSpec spec;
    if (line == Team7BzProbeLine::A2B2)
    {
        spec.y_mm = 144.0;
        spec.z_mm = 34.0;
        spec.csv_basename = "TEAM7_Bz_A2_B2_reference_mT.csv";
    }
    return spec;
}

inline Team7BzProbeLine parseTeam7BzProbeLine(const std::string &name)
{
    if (name == "a2-b2" || name == "a2b2" || name == "A2-B2")
    {
        return Team7BzProbeLine::A2B2;
    }
    return Team7BzProbeLine::A1B1;
}

inline std::string team7BzProbeLineName(Team7BzProbeLine line)
{
    return line == Team7BzProbeLine::A2B2 ? "a2-b2" : "a1-b1";
}

struct Team7BzProbePoint
{
    Real x_mm = 0.0;
    Real bz_ref_phase0_mT = 0.0;
    Real bz_ref_phase90_mT = 0.0;
    Vecd position_m = Vecd::Zero();
    Real bz_sim_mT = 0.0;
};

/** TEAM7 / muFEM Jy probe line slightly below plate top (COMSOL table2 line). */
struct Team7JeyReferenceLineMm
{
    static constexpr Real y_mm = 72.0;
    static constexpr Real z_mm = 18.99;
    static constexpr Real x_start_mm = 0.0;
    static constexpr Real x_end_mm = 288.0;
};

struct Team7JeyProbePoint
{
    Real x_mm = 0.0;
    Real jy_ref_phase0_Am2 = 0.0;
    Real jy_ref_phase90_Am2 = 0.0;
    Real jy_sim_phase0_Am2 = 0.0;
    Real jy_sim_phase90_Am2 = 0.0;
    Vecd position_m = Vecd::Zero();
};

struct Team7BzCompareMetrics
{
    size_t n_probes = 0;
    Real rms_rel_error = 0.0;
    /** RMS of (sim-ref) in mT; diagnostic for phase90 hard gate (future). */
    Real rms_abs_error_mT = 0.0;
    Real max_abs_error_mT = 0.0;
    Real max_rel_error = 0.0;
    size_t peak_ref_index = 0;
    size_t peak_sim_index = 0;
    Real peak_ref_mT = 0.0;
    Real peak_sim_mT = 0.0;
    bool passed = false;
};

struct Team7PlateBzNearProbeSample
{
    Real x_mm = 0.0;
    Real min_dist_mm = 0.0;
    Real plate_bz_mT = 0.0;
    Real host_probe_bz_mT = 0.0;
};

inline std::string resolveTeam7BzReferenceCsvPath(const std::string &reference_dir, Team7BzProbeLine line)
{
    const Team7BzProbeLineSpec spec = team7BzProbeLineSpec(line);
    const StdVec<std::string> candidates = {
        reference_dir + "/" + spec.csv_basename,
        reference_dir + "/team7/" + spec.csv_basename,
        std::string("tests/extra_source_and_tests/3d_examples/reference_data/team7/") + spec.csv_basename,
    };
    for (const std::string &path : candidates)
    {
        if (!path.empty() && fs::exists(path))
        {
            return path;
        }
    }
    return candidates.front();
}

inline std::string resolveTeam7ReferenceCsvPath(const std::string &reference_dir)
{
    return resolveTeam7BzReferenceCsvPath(reference_dir, Team7BzProbeLine::A1B1);
}

inline bool loadTeam7BzReference(const std::string &reference_dir, Team7BzProbeLine line, Real frequency_hz,
                                 Real mm_to_m, StdVec<Team7BzProbePoint> &probes)
{
    probes.clear();
    const Team7BzProbeLineSpec spec = team7BzProbeLineSpec(line);
    const std::string csv_path = resolveTeam7BzReferenceCsvPath(reference_dir, line);
    std::ifstream input(csv_path);
    if (!input)
    {
        std::cout << "[ophelie] TEAM7 Bz reference CSV not found: " << csv_path << std::endl;
        return false;
    }
    const bool use_200hz = frequency_hz > Real(75);
    const size_t phase0_col = use_200hz ? 4 : 2;
    const size_t phase90_col = use_200hz ? 5 : 3;

    std::string line_str;
    std::getline(input, line_str);
    while (std::getline(input, line_str))
    {
        if (line_str.empty() || line_str[0] == '#')
        {
            continue;
        }
        std::stringstream ss(line_str);
        std::string token;
        StdVec<std::string> fields;
        while (std::getline(ss, token, ','))
        {
            fields.push_back(token);
        }
        if (fields.size() <= phase90_col)
        {
            continue;
        }
        Team7BzProbePoint probe;
        probe.x_mm = std::atof(fields[1].c_str());
        probe.bz_ref_phase0_mT = std::atof(fields[phase0_col].c_str());
        probe.bz_ref_phase90_mT = std::atof(fields[phase90_col].c_str());
        probe.position_m = Vecd(probe.x_mm * mm_to_m, spec.y_mm * mm_to_m, spec.z_mm * mm_to_m);
        probes.push_back(probe);
    }
    if (!probes.empty())
    {
        std::cout << "[ophelie] TEAM7 Bz reference loaded: " << csv_path << " line=" << team7BzProbeLineName(line)
                  << " f=" << frequency_hz << " Hz n=" << probes.size() << std::endl;
    }
    return !probes.empty();
}

inline bool loadTeam7BzA1B1Reference(const std::string &reference_dir, Real y_mm, Real z_mm, Real mm_to_m,
                                     StdVec<Team7BzProbePoint> &probes)
{
    (void)y_mm;
    (void)z_mm;
    return loadTeam7BzReference(reference_dir, Team7BzProbeLine::A1B1, Real(50), mm_to_m, probes);
}

inline std::string resolveTeam7JeyReferenceCsvPath(const std::string &reference_dir)
{
    const StdVec<std::string> candidates = {
        reference_dir + "/TEAM7_Jey_A1_surface_reference_Am2.csv",
        reference_dir + "/team7/TEAM7_Jey_A1_surface_reference_Am2.csv",
        "tests/extra_source_and_tests/3d_examples/reference_data/team7/TEAM7_Jey_A1_surface_reference_Am2.csv",
    };
    for (const std::string &path : candidates)
    {
        if (!path.empty() && fs::exists(path))
        {
            return path;
        }
    }
    return candidates.front();
}

/** Optional COMSOL/muFEM Jy line reference; returns false if CSV not bundled yet. */
inline bool loadTeam7JeyReference(const std::string &reference_dir, Real frequency_hz, Real mm_to_m,
                                  StdVec<Team7JeyProbePoint> &probes)
{
    probes.clear();
    const std::string csv_path = resolveTeam7JeyReferenceCsvPath(reference_dir);
    std::ifstream input(csv_path);
    if (!input)
    {
        std::cout << "[ophelie] TEAM7 Jey reference CSV not found (diagnostic sim-only): " << csv_path << std::endl;
        return false;
    }
    const bool use_200hz = frequency_hz > Real(75);
    const size_t phase0_col = use_200hz ? 4 : 2;
    const size_t phase90_col = use_200hz ? 5 : 3;

    std::string line_str;
    std::getline(input, line_str);
    while (std::getline(input, line_str))
    {
        if (line_str.empty() || line_str[0] == '#')
        {
            continue;
        }
        std::stringstream ss(line_str);
        std::string token;
        StdVec<std::string> fields;
        while (std::getline(ss, token, ','))
        {
            fields.push_back(token);
        }
        if (fields.size() <= phase90_col)
        {
            continue;
        }
        Team7JeyProbePoint probe;
        probe.x_mm = std::atof(fields[1].c_str());
        probe.jy_ref_phase0_Am2 = std::atof(fields[phase0_col].c_str());
        probe.jy_ref_phase90_Am2 = std::atof(fields[phase90_col].c_str());
        probe.position_m =
            Vecd(probe.x_mm * mm_to_m, Team7JeyReferenceLineMm::y_mm * mm_to_m, Team7JeyReferenceLineMm::z_mm * mm_to_m);
        probes.push_back(probe);
    }
    if (!probes.empty())
    {
        std::cout << "[ophelie] TEAM7 Jey reference loaded: " << csv_path << " f=" << frequency_hz
                  << " Hz n=" << probes.size() << std::endl;
    }
    return !probes.empty();
}

inline bool team7JeyReferenceHasFrequencyData(const StdVec<Team7JeyProbePoint> &probes, Real frequency_hz)
{
    const bool use_200hz = frequency_hz > Real(75);
    for (const Team7JeyProbePoint &probe : probes)
    {
        if (use_200hz)
        {
            if (std::abs(probe.jy_ref_phase0_Am2) > TinyReal || std::abs(probe.jy_ref_phase90_Am2) > TinyReal)
            {
                return true;
            }
        }
        else if (std::abs(probe.jy_ref_phase0_Am2) > TinyReal || std::abs(probe.jy_ref_phase90_Am2) > TinyReal)
        {
            return true;
        }
    }
    return false;
}

inline void buildTeam7JeyProbeLineFromBzReference(const StdVec<Team7BzProbePoint> &bz_reference_probes, Real mm_to_m,
                                                  StdVec<Team7JeyProbePoint> &jey_probes)
{
    jey_probes.clear();
    for (const Team7BzProbePoint &bz_probe : bz_reference_probes)
    {
        Team7JeyProbePoint jey_probe;
        jey_probe.x_mm = bz_probe.x_mm;
        jey_probe.position_m =
            Vecd(bz_probe.x_mm * mm_to_m, Team7JeyReferenceLineMm::y_mm * mm_to_m, Team7JeyReferenceLineMm::z_mm * mm_to_m);
        jey_probes.push_back(jey_probe);
    }
}

/** Host Biot–Savart: coil J at probe points in air (TEAM7 A1–B1 line, not on PlateBody). */
inline void evaluateCoilBiotSavartBzAtProbes(BaseParticles &coil_particles, const OphelieCoilFieldNames &coil_names,
                                             const OphelieParameters &params, const StdVec<Team7BzProbePoint> &probes,
                                             StdVec<Team7BzProbePoint> &results)
{
    results = probes;
    syncVariableToHost<Vecd>(coil_particles, "Position");
    syncVariableToHost<Vecd>(coil_particles, coil_names.j_src_real);
    syncVariableToHost<Real>(coil_particles, "VolumetricMeasure");

    const Vecd *coil_pos = coil_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j_src = coil_particles.getVariableDataByName<Vecd>(coil_names.j_src_real);
    const Real *vol = coil_particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n_coil = coil_particles.TotalRealParticles();
    const Real coeff = params.mu0_ / (4.0 * Pi);
    const Real eps2 = params.softening_length_ * params.softening_length_;

    for (Team7BzProbePoint &probe : results)
    {
        Vecd b_sum = Vecd::Zero();
        for (size_t j = 0; j < n_coil; ++j)
        {
            const Vecd r = probe.position_m - coil_pos[j];
            const Real r2 = r.squaredNorm() + eps2;
            const Real inv_r = 1.0 / std::sqrt(r2);
            const Real inv_r3 = inv_r / r2;
            const Vecd jv = j_src[j] * vol[j];
            b_sum += coeff * jv.cross(r) * inv_r3;
        }
        probe.bz_sim_mT = b_sum[2] * 1000.0;
    }
}

inline Team7BzCompareMetrics compareTeam7ScalarAgainstReference(const StdVec<Real> &sim_values,
                                                                const StdVec<Real> &ref_values,
                                                                Real smoke_rms_threshold = 0.5)
{
    Team7BzCompareMetrics metrics;
    metrics.n_probes = std::min(sim_values.size(), ref_values.size());
    Real sum_ref_sq = 0.0;
    Real sum_err_sq = 0.0;
    for (size_t i = 0; i < metrics.n_probes; ++i)
    {
        const Real ref = ref_values[i];
        const Real sim = sim_values[i];
        const Real err = sim - ref;
        sum_ref_sq += ref * ref;
        sum_err_sq += err * err;
        metrics.max_abs_error_mT = std::max(metrics.max_abs_error_mT, std::abs(err));
        metrics.max_rel_error = std::max(metrics.max_rel_error, std::abs(err) / (std::abs(ref) + TinyReal));
        if (std::abs(ref) >= std::abs(metrics.peak_ref_mT))
        {
            metrics.peak_ref_mT = ref;
            metrics.peak_ref_index = i;
        }
        if (std::abs(sim) >= std::abs(metrics.peak_sim_mT))
        {
            metrics.peak_sim_mT = sim;
            metrics.peak_sim_index = i;
        }
    }
    metrics.rms_abs_error_mT =
        metrics.n_probes > 0 ? std::sqrt(sum_err_sq / static_cast<Real>(metrics.n_probes)) : Real(0);
    metrics.rms_rel_error = std::sqrt(sum_err_sq) / (std::sqrt(sum_ref_sq) + TinyReal);
    metrics.passed = metrics.n_probes > 0 && metrics.rms_rel_error < smoke_rms_threshold;
    return metrics;
}

inline Team7BzCompareMetrics compareTeam7BzAgainstReference(const StdVec<Team7BzProbePoint> &probes,
                                                            const StdVec<Real> &sim_mT,
                                                            const StdVec<Real> &ref_mT, Real smoke_rms_threshold = 0.5)
{
    Team7BzCompareMetrics metrics;
    metrics.n_probes = std::min({probes.size(), sim_mT.size(), ref_mT.size()});
    Real sum_ref_sq = 0.0;
    Real sum_err_sq = 0.0;
    for (size_t i = 0; i < metrics.n_probes; ++i)
    {
        const Real ref = ref_mT[i];
        const Real sim = sim_mT[i];
        const Real err = sim - ref;
        sum_ref_sq += ref * ref;
        sum_err_sq += err * err;
        metrics.max_abs_error_mT = std::max(metrics.max_abs_error_mT, std::abs(err));
        metrics.max_rel_error = std::max(metrics.max_rel_error, std::abs(err) / (std::abs(ref) + TinyReal));
        if (std::abs(ref) >= std::abs(metrics.peak_ref_mT))
        {
            metrics.peak_ref_mT = ref;
            metrics.peak_ref_index = i;
        }
        if (std::abs(sim) >= std::abs(metrics.peak_sim_mT))
        {
            metrics.peak_sim_mT = sim;
            metrics.peak_sim_index = i;
        }
    }
    metrics.rms_abs_error_mT =
        metrics.n_probes > 0 ? std::sqrt(sum_err_sq / static_cast<Real>(metrics.n_probes)) : Real(0);
    metrics.rms_rel_error = std::sqrt(sum_err_sq) / (std::sqrt(sum_ref_sq) + TinyReal);
    metrics.passed = metrics.n_probes > 0 && metrics.rms_rel_error < smoke_rms_threshold;
    return metrics;
}

inline Team7BzCompareMetrics compareTeam7BzPhase0(const StdVec<Team7BzProbePoint> &probes, Real smoke_rms_threshold = 0.5)
{
    Team7BzCompareMetrics metrics;
    metrics.n_probes = probes.size();
    Real sum_ref_sq = 0.0;
    Real sum_err_sq = 0.0;
    for (size_t i = 0; i < probes.size(); ++i)
    {
        const Real ref = probes[i].bz_ref_phase0_mT;
        const Real sim = probes[i].bz_sim_mT;
        const Real err = sim - ref;
        sum_ref_sq += ref * ref;
        sum_err_sq += err * err;
        metrics.max_abs_error_mT = std::max(metrics.max_abs_error_mT, std::abs(err));
        metrics.max_rel_error = std::max(metrics.max_rel_error, std::abs(err) / (std::abs(ref) + TinyReal));
        if (std::abs(ref) >= std::abs(metrics.peak_ref_mT))
        {
            metrics.peak_ref_mT = ref;
            metrics.peak_ref_index = i;
        }
        if (std::abs(sim) >= std::abs(metrics.peak_sim_mT))
        {
            metrics.peak_sim_mT = sim;
            metrics.peak_sim_index = i;
        }
    }
    metrics.rms_abs_error_mT =
        metrics.n_probes > 0 ? std::sqrt(sum_err_sq / static_cast<Real>(metrics.n_probes)) : Real(0);
    metrics.rms_rel_error = std::sqrt(sum_err_sq) / (std::sqrt(sum_ref_sq) + TinyReal);
    metrics.passed = metrics.n_probes > 0 && metrics.rms_rel_error < smoke_rms_threshold;
    return metrics;
}

/** Nearest PlateBody particle to each air probe (3D); sanity-check device Biot vs host probe. */
inline void sampleNearestPlateCoilBzAtProbes(BaseParticles &plate_particles, const OphelieGlassFieldNames &plate_names,
                                              const StdVec<Team7BzProbePoint> &host_probes, Real mm_to_m,
                                              StdVec<Team7PlateBzNearProbeSample> &samples)
{
    samples.clear();
    samples.resize(host_probes.size());
    syncVariableToHost<Vecd>(plate_particles, "Position");
    syncVariableToHost<Vecd>(plate_particles, plate_names.b_coil_real);
    const Vecd *pos = plate_particles.getVariableDataByName<Vecd>("Position");
    const Vecd *b_coil = plate_particles.getVariableDataByName<Vecd>(plate_names.b_coil_real);
    const size_t n_plate = plate_particles.TotalRealParticles();

    for (size_t ip = 0; ip < host_probes.size(); ++ip)
    {
        samples[ip].x_mm = host_probes[ip].x_mm;
        samples[ip].host_probe_bz_mT = host_probes[ip].bz_sim_mT;
        Real min_dist_m = std::numeric_limits<Real>::infinity();
        Real plate_bz_mT = 0.0;
        const Vecd &probe = host_probes[ip].position_m;
        for (size_t i = 0; i < n_plate; ++i)
        {
            const Real dist = (pos[i] - probe).norm();
            if (dist < min_dist_m)
            {
                min_dist_m = dist;
                plate_bz_mT = b_coil[i][2] * 1000.0;
            }
        }
        samples[ip].min_dist_mm = min_dist_m / mm_to_m;
        samples[ip].plate_bz_mT = plate_bz_mT;
    }
}

inline void writeTeam7JeyProbeCsv(const std::string &output_path, const StdVec<Team7JeyProbePoint> &probes,
                                bool reference_ok)
{
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 Jey probe CSV: " << output_path << std::endl;
        return;
    }
    if (reference_ok)
    {
        out << "x_mm,Jy_ref_phase0_Am2,Jy_ref_phase90_Am2,Jy_sim_phase0_Am2,Jy_sim_phase90_Am2\n";
        for (const Team7JeyProbePoint &probe : probes)
        {
            out << probe.x_mm << "," << probe.jy_ref_phase0_Am2 << "," << probe.jy_ref_phase90_Am2 << ","
                << probe.jy_sim_phase0_Am2 << "," << probe.jy_sim_phase90_Am2 << "\n";
        }
    }
    else
    {
        out << "x_mm,Jy_sim_phase0_Am2,Jy_sim_phase90_Am2\n";
        for (const Team7JeyProbePoint &probe : probes)
        {
            out << probe.x_mm << "," << probe.jy_sim_phase0_Am2 << "," << probe.jy_sim_phase90_Am2 << "\n";
        }
    }
    std::cout << "[ophelie] TEAM7 Jey probe CSV: " << output_path << std::endl;
}

inline void writeTeam7BzCompareCsv(const std::string &output_path, const StdVec<Team7BzProbePoint> &probes)
{
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    if (!out)
    {
        std::cout << "[ophelie] could not write TEAM7 Bz compare CSV: " << output_path << std::endl;
        return;
    }
    out << "x_mm,Bz_ref_phase0_mT,Bz_sim_mT,err_mT,rel_err\n";
    for (const Team7BzProbePoint &probe : probes)
    {
        const Real err = probe.bz_sim_mT - probe.bz_ref_phase0_mT;
        const Real rel = err / (std::abs(probe.bz_ref_phase0_mT) + TinyReal);
        out << probe.x_mm << "," << probe.bz_ref_phase0_mT << "," << probe.bz_sim_mT << "," << err << "," << rel << "\n";
    }
    std::cout << "[ophelie] TEAM7 Bz compare CSV: " << output_path << std::endl;
}

/** Circular-filament diagnostic: N*I on a loop at coil mean radius (not TEAM7 racetrack). */
inline void evaluateTeam7CircularLoopBiotSavartBzAtProbes(const Vecd &loop_center_m, Real loop_radius_m,
                                                          Real total_current_a, Real mu0,
                                                          const StdVec<Team7BzProbePoint> &probes,
                                                          StdVec<Team7BzProbePoint> &results, size_t n_segments = 720)
{
    results = probes;
    const Real coeff = mu0 / (4.0 * Pi);
    for (Team7BzProbePoint &probe : results)
    {
        Vecd b_sum = Vecd::Zero();
        for (size_t k = 0; k < n_segments; ++k)
        {
            const Real phi0 = 2.0 * Pi * static_cast<Real>(k) / static_cast<Real>(n_segments);
            const Real phi1 = 2.0 * Pi * static_cast<Real>(k + 1) / static_cast<Real>(n_segments);
            const Vecd p0(loop_center_m[0] + loop_radius_m * std::cos(phi0),
                          loop_center_m[1] + loop_radius_m * std::sin(phi0), loop_center_m[2]);
            const Vecd p1(loop_center_m[0] + loop_radius_m * std::cos(phi1),
                          loop_center_m[1] + loop_radius_m * std::sin(phi1), loop_center_m[2]);
            const Vecd dl = p1 - p0;
            const Vecd r = probe.position_m - p0;
            const Real r2 = r.squaredNorm() + TinyReal;
            const Real inv_r3 = 1.0 / (std::sqrt(r2) * r2);
            b_sum += coeff * total_current_a * dl.cross(r) * inv_r3;
        }
        probe.bz_sim_mT = b_sum[2] * 1000.0;
    }
}

/** Rectangular filament at coil STL xy footprint (racetrack diagnostic, z = loop_z). */
inline void evaluateTeam7RectangularLoopBiotSavartBzAtProbes(const Team7NativeStlMeshBBoxMm &mesh_bbox, Real loop_z_mm,
                                                             Real total_current_a, Real mu0,
                                                             const StdVec<Team7BzProbePoint> &probes,
                                                             StdVec<Team7BzProbePoint> &results,
                                                             size_t points_per_edge = 160)
{
    results = probes;
    const Real coeff = mu0 / (4.0 * Pi);
    const Real mm_to_m = 1.0e-3;
    const Real x0 = mesh_bbox.coil_lower_[0];
    const Real x1 = mesh_bbox.coil_upper_[0];
    const Real y0 = mesh_bbox.coil_lower_[1];
    const Real y1 = mesh_bbox.coil_upper_[1];
    const Real z_m = loop_z_mm * mm_to_m;
    const Vec3d corners[4] = {{x0, y0, loop_z_mm}, {x1, y0, loop_z_mm}, {x1, y1, loop_z_mm}, {x0, y1, loop_z_mm}};
    StdVec<Vecd> polyline;
    polyline.reserve(4 * points_per_edge);
    for (int edge = 0; edge < 4; ++edge)
    {
        const Vec3d a = corners[edge];
        const Vec3d b = corners[(edge + 1) % 4];
        for (size_t k = 0; k < points_per_edge; ++k)
        {
            const Real t = static_cast<Real>(k) / static_cast<Real>(points_per_edge);
            polyline.push_back(Vecd((a[0] + t * (b[0] - a[0])) * mm_to_m, (a[1] + t * (b[1] - a[1])) * mm_to_m, z_m));
        }
    }
    const size_t n_pts = polyline.size();
    for (Team7BzProbePoint &probe : results)
    {
        Vecd b_sum = Vecd::Zero();
        for (size_t k = 0; k < n_pts; ++k)
        {
            const Vecd p0 = polyline[k];
            const Vecd p1 = polyline[(k + 1) % n_pts];
            const Vecd dl = p1 - p0;
            const Vecd r = probe.position_m - p0;
            const Real r2 = r.squaredNorm() + TinyReal;
            const Real inv_r3 = 1.0 / (std::sqrt(r2) * r2);
            b_sum += coeff * total_current_a * dl.cross(r) * inv_r3;
        }
        probe.bz_sim_mT = b_sum[2] * 1000.0;
    }
}

inline Team7BzCompareMetrics compareTeam7BzPhase0OverXRange(const StdVec<Team7BzProbePoint> &probes, Real x_min_mm,
                                                            Real x_max_mm, Real smoke_rms_threshold = 0.5)
{
    StdVec<Team7BzProbePoint> subset;
    for (const Team7BzProbePoint &probe : probes)
    {
        if (probe.x_mm >= x_min_mm && probe.x_mm <= x_max_mm)
        {
            subset.push_back(probe);
        }
    }
    return compareTeam7BzPhase0(subset, smoke_rms_threshold);
}

inline void printTeam7BzCompareReport(const StdVec<Team7BzProbePoint> &probes, const Team7BzCompareMetrics &metrics,
                                      Real smoke_rms_threshold)
{
    std::cout << "[ophelie] TEAM7 Bz A1-B1 compare n=" << metrics.n_probes << " rms_rel_err=" << metrics.rms_rel_error
              << " max_abs_err_mT=" << metrics.max_abs_error_mT << " max_rel_err=" << metrics.max_rel_error
              << " threshold=" << smoke_rms_threshold << " passed=" << (metrics.passed ? 1 : 0) << std::endl;
    if (metrics.n_probes > 0)
    {
        std::cout << "[ophelie]   peak_ref x_mm=" << probes[metrics.peak_ref_index].x_mm
                  << " |Bz_ref|_mT=" << std::abs(metrics.peak_ref_mT) << " peak_sim x_mm="
                  << probes[metrics.peak_sim_index].x_mm << " |Bz_sim|_mT=" << std::abs(metrics.peak_sim_mT)
                  << std::endl;
    }
    std::cout << "[ophelie]   x_mm  Bz_ref_mT  Bz_sim_mT  err_mT  rel_err" << std::endl;
    for (size_t i = 0; i < probes.size(); ++i)
    {
        const Real err = probes[i].bz_sim_mT - probes[i].bz_ref_phase0_mT;
        const Real rel = err / (std::abs(probes[i].bz_ref_phase0_mT) + TinyReal);
        std::cout << "[ophelie]   " << probes[i].x_mm << "  " << probes[i].bz_ref_phase0_mT << "  "
                  << probes[i].bz_sim_mT << "  " << err << "  " << rel << std::endl;
    }
}

inline void printTeam7PlateBzNearProbeReport(const StdVec<Team7PlateBzNearProbeSample> &samples)
{
    std::cout << "[ophelie] nearest PlateBody Bz vs air probe (expect large min_dist_mm; values are diagnostic only)"
              << std::endl;
    for (const Team7PlateBzNearProbeSample &sample : samples)
    {
        std::cout << "[ophelie]   x_mm=" << sample.x_mm << " min_dist_mm=" << sample.min_dist_mm
                  << " plate_Bz_mT=" << sample.plate_bz_mT << " host_probe_Bz_mT=" << sample.host_probe_bz_mT
                  << std::endl;
    }
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_PROBE_H
