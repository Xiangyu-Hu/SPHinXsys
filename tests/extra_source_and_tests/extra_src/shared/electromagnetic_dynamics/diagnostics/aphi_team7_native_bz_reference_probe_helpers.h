#ifndef APHI_TEAM7_NATIVE_BZ_REFERENCE_PROBE_HELPERS_H
#define APHI_TEAM7_NATIVE_BZ_REFERENCE_PROBE_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_reload_geometry_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_solve_field_audit_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_mr_probe_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_probe_reference_alignment_helpers.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** muFEM/COMSOL TEAM7 probe lines in SI [m] (official TEAM7 diagram values × 0.001). */
struct AphiTeam7NativeProbeLineDefinition
{
    std::string probe_id;
    Real y_m = 0.072;
    Real z_m = 0.034;
    Real x_start_m = 0.0;
    Real x_end_m = 0.288;
    Real x_step_m = 0.018;
    std::string reference_csv_path;
};

struct AphiTeam7BzReferenceRow
{
    Real x_mm = 0.0;
    Real bz_50hz_phase0_mT = 0.0;
    Real bz_50hz_phase90_mT = 0.0;
};

struct AphiTeam7BzSignConventionAudit
{
    int sign_bz_real = -1;
    int sign_bz_imag = 1;
    Real tesla_to_mT_scale = 1000.0;
    Real profile_rel_err_real = 0.0;
    Real profile_rel_err_imag = 0.0;
    Real profile_rel_err_sum = 0.0;
    Real max_point_rel_err_real = 0.0;
    Real max_point_rel_err_imag = 0.0;
    Real worst_point_rel_err_x_mm = 0.0;
};

struct AphiTeam7NativeBzProbeLineResult
{
    AphiTeam7NativeProbeLineDefinition definition{};
    StdVec<AphiTeam7BzReferenceRow> reference{};
    StdVec<AphiProbeLineSample> samples{};
    /** Active sign convention for CSV / pass gate (frozen or best-scan). */
    AphiTeam7BzSignConventionAudit sign_audit{};
    /** Best-of-four scan; filled when frozen sign mode runs diagnostic comparison. */
    AphiTeam7BzSignConventionAudit sign_scan_audit{};
    bool sign_scan_audit_valid = false;
    Real max_nearest_distance_mm = 0.0;
    bool reference_loaded = false;
    bool probe_csv_written = false;
    bool comparison_csv_written = false;
};

struct AphiTeam7NativeBzReferenceProbeSpec
{
    AphiTeam7NativeReloadSourceDrivenSmokeSpec solve_spec{};
    std::string report_dir = "team7_native_bz_reference_report";
    /** Max probe-to-particle distance [m] (default 3× dp_0 from reload metadata). */
    Real max_probe_nearest_distance_m = 3.0 * AphiTeam7NativeGeometryConstants::dp_0_default;
    /** P3: profile L2 gate on A1-B1 @ 50 Hz (under frozen TEAM7 sign). */
    Real max_a1b1_profile_rel_err = 0.25;
    Real tesla_to_mT_scale = 1000.0;
    /** P2/P3: require |Bz| on A1-B1 in [min, max] mT before profile error gates tighten. */
    Real min_max_abs_bz_a1b1_mT = 0.01;
    Real max_max_abs_bz_a1b1_mT = 100.0;
    bool require_bz_magnitude_sanity = true;
    /** P2 vacuum: compare mT magnitude only, not profile vs eddy-current reference. */
    bool require_a1b1_profile_gate = true;
    /**
     * P4: freeze muFEM-like convention for conductive TEAM7 (scan on A1-B1 gave (+real,+imag) vs ref).
     * Maps our CK Bz to mT: phase0_mT = sign_bz_real * Bz_real[T]*1e3, phase90_mT = sign_bz_imag * Bz_imag[T]*1e3.
     */
    bool use_frozen_team7_bz_sign = false;
    int frozen_sign_bz_real = 1;
    int frozen_sign_bz_imag = -1;
    bool run_sign_scan_diagnostic = true;
    bool require_frozen_sign_matches_scan = true;
    /** TEAM7 A1-B1/A2-B2 lie in air; default sample nearest air particle only. */
    bool probe_b_sample_air_particles_only = true;
    /** After curl-B probe: write Air/Coil/Plate VTP under ./output/ (CLI: --team7-write-vtp). */
    bool write_vtp = false;
};

struct AphiTeam7NativeBzReferenceProbeSummary
{
    bool reload_files_present = false;
    bool reference_data_present = false;
    bool solver_ok = false;
    bool vtp_written = false;
    Real final_true_rel = 0.0;
    AphiTeam7NativeBzProbeLineResult a1_b1{};
    AphiTeam7NativeBzProbeLineResult a2_b2{};
    bool summary_csv_written = false;
};

inline StdVec<Vecd> buildTeam7NativeProbeLinePositionsSi(const AphiTeam7NativeProbeLineDefinition &definition)
{
    StdVec<Vecd> positions;
    for (Real x_m = definition.x_start_m; x_m <= definition.x_end_m + TinyReal; x_m += definition.x_step_m)
    {
        positions.push_back(Vecd(x_m, definition.y_m, definition.z_m));
    }
    return positions;
}

inline bool loadTeam7BzReferenceCsv(const std::string &path, StdVec<AphiTeam7BzReferenceRow> &rows)
{
    std::ifstream input(path);
    if (!input)
    {
        return false;
    }
    std::string header;
    if (!std::getline(input, header))
    {
        return false;
    }
    rows.clear();
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty())
        {
            continue;
        }
        std::stringstream ss(line);
        std::string token;
        StdVec<std::string> fields;
        while (std::getline(ss, token, ','))
        {
            fields.push_back(token);
        }
        if (fields.size() < 4)
        {
            continue;
        }
        AphiTeam7BzReferenceRow row;
        row.x_mm = std::stod(fields[1]);
        row.bz_50hz_phase0_mT = std::stod(fields[2]);
        row.bz_50hz_phase90_mT = std::stod(fields[3]);
        rows.push_back(row);
    }
    return !rows.empty();
}

inline void execCurlBOnNativeReloadCase(AphiTeam7NativeReloadContactCase &case_setup, const AphiVariableNames &names)
{
    static constexpr const char *k_b_real = "NativeProbeBReal";
    static constexpr const char *k_b_imag = "NativeProbeBImag";
    execBodyCurlBFromAdaptiveInnerDiagnostic<Team7NativeAdaptiveAirInner, Team7NativeAdaptiveAirBody>(
        case_setup.air_body(), case_setup.air_inner(), names, k_b_real, k_b_imag, AphiBCurlDiagnosticMode::BCorrectedGrad);
    execBodyCurlBFromAdaptiveInnerDiagnostic<Team7NativeAdaptiveCoilInner, Team7NativeAdaptiveCoilBody>(
        case_setup.coil_body(), case_setup.coil_inner(), names, k_b_real, k_b_imag, AphiBCurlDiagnosticMode::BCorrectedGrad);
    execBodyCurlBFromAdaptiveInnerDiagnostic<Team7NativeAdaptivePlateInner, Team7NativeAdaptivePlateBody>(
        case_setup.plate_body(), case_setup.plate_inner(), names, k_b_real, k_b_imag,
        AphiBCurlDiagnosticMode::BCorrectedGrad);
}

inline StdVec<AphiProbeLineSample> sampleBProbeFromNativeReloadCase(AphiTeam7NativeReloadContactCase &case_setup,
                                                                    const StdVec<Vecd> &sample_positions,
                                                                    bool air_particles_only = true)
{
    static constexpr const char *k_b_real = "NativeProbeBReal";
    static constexpr const char *k_b_imag = "NativeProbeBImag";
    BaseParticles &air_particles = case_setup.air_body().getBaseParticles();
    const Kernel &kernel = *case_setup.air_body().getSPHAdaptation().getKernel();
    const AphiTeam7NativeGeometryConfig &geometry_config = case_setup.geometry_config;

    StdVec<AphiProbeLineSample> samples;
    samples.reserve(sample_positions.size());
    if (air_particles_only)
    {
        for (const Vecd &sample_position : sample_positions)
        {
            samples.push_back(sampleNearestAirBWithTeam7KernelSupport(air_particles, kernel, geometry_config,
                                                                      sample_position, k_b_real, k_b_imag));
        }
        return samples;
    }

    BaseParticles *coil_particles = &case_setup.coil_body().getBaseParticles();
    BaseParticles *plate_particles = &case_setup.plate_body().getBaseParticles();
    BaseParticles *const all_views[] = {&air_particles, coil_particles, plate_particles};
    for (const Vecd &sample_position : sample_positions)
    {
        AphiProbeLineSample sample;
        sample.position = sample_position;
        Real best_distance_squared = std::numeric_limits<Real>::max();
        for (size_t v = 0; v != 3; ++v)
        {
            BaseParticles *particles = all_views[v];
            syncVariableToHost<Vecd>(*particles, "Position");
            syncVariableToHost<Vecd>(*particles, k_b_real);
            syncVariableToHost<Vecd>(*particles, k_b_imag);
            const Vecd *positions = particles->getVariableDataByName<Vecd>("Position");
            const Vecd *b_real = particles->getVariableDataByName<Vecd>(k_b_real);
            const Vecd *b_imag = particles->getVariableDataByName<Vecd>(k_b_imag);
            const size_t total_real_particles = particles->TotalRealParticles();
            for (size_t i = 0; i != total_real_particles; ++i)
            {
                const Real distance_squared = (positions[i] - sample_position).squaredNorm();
                if (distance_squared < best_distance_squared)
                {
                    best_distance_squared = distance_squared;
                    sample.nearest_particle_index = i;
                    sample.neighbor_distance = std::sqrt(distance_squared);
                    sample.b_real = b_real[i];
                    sample.b_imag = b_imag[i];
                    sample.b_abs = vec3Abs(sample.b_real);
                }
            }
        }
        samples.push_back(sample);
    }
    return samples;
}

inline Real profileL2RelativeDifference(const StdVec<Real> &reference, const StdVec<Real> &approximate)
{
    if (reference.size() != approximate.size() || reference.empty())
    {
        return std::numeric_limits<Real>::infinity();
    }
    Real sum_sq = 0.0;
    Real ref_sum_sq = 0.0;
    for (size_t i = 0; i != reference.size(); ++i)
    {
        const Real diff = approximate[i] - reference[i];
        sum_sq += diff * diff;
        ref_sum_sq += reference[i] * reference[i];
    }
    return std::sqrt(sum_sq) / (std::sqrt(ref_sum_sq) + TinyReal);
}

inline AphiTeam7BzSignConventionAudit buildBzSignConventionAudit(const StdVec<AphiProbeLineSample> &samples,
                                                                 const StdVec<AphiTeam7BzReferenceRow> &reference_rows,
                                                                 Real tesla_to_mT_scale, int sign_bz_real,
                                                                 int sign_bz_imag)
{
    AphiTeam7BzSignConventionAudit audit;
    audit.sign_bz_real = sign_bz_real;
    audit.sign_bz_imag = sign_bz_imag;
    audit.tesla_to_mT_scale = tesla_to_mT_scale;
    StdVec<Real> our_real;
    StdVec<Real> our_imag;
    StdVec<Real> ref_real;
    StdVec<Real> ref_imag;
    const size_t count = std::min(samples.size(), reference_rows.size());
    our_real.reserve(count);
    our_imag.reserve(count);
    ref_real.reserve(count);
    ref_imag.reserve(count);
    for (size_t i = 0; i != count; ++i)
    {
        our_real.push_back(sign_bz_real * tesla_to_mT_scale * samples[i].b_real[2]);
        our_imag.push_back(sign_bz_imag * tesla_to_mT_scale * samples[i].b_imag[2]);
        ref_real.push_back(reference_rows[i].bz_50hz_phase0_mT);
        ref_imag.push_back(reference_rows[i].bz_50hz_phase90_mT);
    }
    audit.profile_rel_err_real = profileL2RelativeDifference(ref_real, our_real);
    audit.profile_rel_err_imag = profileL2RelativeDifference(ref_imag, our_imag);
    audit.profile_rel_err_sum = audit.profile_rel_err_real + audit.profile_rel_err_imag;
    /** Skip near-zero reference samples (e.g. Bz zero-crossing at x≈90 mm) for pointwise rel err. */
    constexpr Real min_ref_mT_for_pointwise_rel = 0.1;
    for (size_t i = 0; i != count; ++i)
    {
        if (std::fabs(ref_real[i]) > min_ref_mT_for_pointwise_rel)
        {
            const Real rel_real =
                std::fabs(our_real[i] - ref_real[i]) / (std::fabs(ref_real[i]) + TinyReal);
            if (rel_real > audit.max_point_rel_err_real)
            {
                audit.max_point_rel_err_real = rel_real;
                audit.worst_point_rel_err_x_mm = reference_rows[i].x_mm;
            }
        }
        if (std::fabs(ref_imag[i]) > min_ref_mT_for_pointwise_rel)
        {
            const Real rel_imag =
                std::fabs(our_imag[i] - ref_imag[i]) / (std::fabs(ref_imag[i]) + TinyReal);
            if (rel_imag > audit.max_point_rel_err_imag)
            {
                audit.max_point_rel_err_imag = rel_imag;
            }
        }
    }
    return audit;
}

inline AphiTeam7BzSignConventionAudit pickBestBzSignConvention(const StdVec<AphiProbeLineSample> &samples,
                                                              const StdVec<AphiTeam7BzReferenceRow> &reference_rows,
                                                              Real tesla_to_mT_scale)
{
    AphiTeam7BzSignConventionAudit best;
    best.profile_rel_err_sum = std::numeric_limits<Real>::infinity();
    const StdVec<int> sign_options = {1, -1};
    for (const int sign_real : sign_options)
    {
        for (const int sign_imag : sign_options)
        {
            const AphiTeam7BzSignConventionAudit candidate =
                buildBzSignConventionAudit(samples, reference_rows, tesla_to_mT_scale, sign_real, sign_imag);
            if (candidate.profile_rel_err_sum < best.profile_rel_err_sum)
            {
                best = candidate;
            }
        }
    }
    return best;
}

inline bool frozenTeam7BzSignMatchesScan(const AphiTeam7BzSignConventionAudit &frozen,
                                         const AphiTeam7BzSignConventionAudit &scanned)
{
    return frozen.sign_bz_real == scanned.sign_bz_real && frozen.sign_bz_imag == scanned.sign_bz_imag;
}

inline bool writeNativeBzProbeCsv(const std::string &path, const StdVec<AphiProbeLineSample> &samples,
                                const AphiTeam7BzSignConventionAudit &sign_audit)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "x_mm,y_mm,z_mm,nearest_distance_mm,Bz_real_T,Bz_imag_T,"
              "our_Bz_real_mT,our_Bz_imag_mT,our_minus_Bz_real_mT,our_minus_Bz_imag_mT\n";
    for (const AphiProbeLineSample &sample : samples)
    {
        const Real our_real_mT = sign_audit.sign_bz_real * sign_audit.tesla_to_mT_scale * sample.b_real[2];
        const Real our_imag_mT = sign_audit.sign_bz_imag * sign_audit.tesla_to_mT_scale * sample.b_imag[2];
        output << (1.0e3 * sample.position[0]) << "," << (1.0e3 * sample.position[1]) << ","
               << (1.0e3 * sample.position[2]) << "," << (1.0e3 * sample.neighbor_distance) << "," << sample.b_real[2]
               << "," << sample.b_imag[2] << "," << our_real_mT
               << "," << our_imag_mT << "," << (-our_real_mT) << "," << (-our_imag_mT) << "\n";
    }
    return true;
}

inline bool writeNativeBzComparisonCsv(const std::string &path, const StdVec<AphiProbeLineSample> &samples,
                                       const StdVec<AphiTeam7BzReferenceRow> &reference_rows,
                                       const AphiTeam7BzSignConventionAudit &sign_audit)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "x_mm,our_Bz_real_mT,our_Bz_imag_mT,our_minus_Bz_real_mT,our_minus_Bz_imag_mT,"
              "ref_Bz_50Hz_phase0_mT,ref_Bz_50Hz_phase90_mT,nearest_distance_mm\n";
    const size_t count = std::min(samples.size(), reference_rows.size());
    for (size_t i = 0; i != count; ++i)
    {
        const Real our_real_mT = sign_audit.sign_bz_real * sign_audit.tesla_to_mT_scale * samples[i].b_real[2];
        const Real our_imag_mT = sign_audit.sign_bz_imag * sign_audit.tesla_to_mT_scale * samples[i].b_imag[2];
        output << reference_rows[i].x_mm << "," << our_real_mT << "," << our_imag_mT << "," << (-our_real_mT) << ","
               << (-our_imag_mT) << "," << reference_rows[i].bz_50hz_phase0_mT << ","
               << reference_rows[i].bz_50hz_phase90_mT << "," << (1.0e3 * samples[i].neighbor_distance) << "\n";
    }
    return true;
}

inline AphiTeam7NativeBzProbeLineResult runNativeBzProbeLine(AphiTeam7NativeReloadContactCase &case_setup,
                                                             const AphiTeam7NativeProbeLineDefinition &definition,
                                                             const AphiTeam7NativeBzReferenceProbeSpec &spec)
{
    AphiTeam7NativeBzProbeLineResult result;
    result.definition = definition;
    result.reference_loaded = loadTeam7BzReferenceCsv(definition.reference_csv_path, result.reference);
    if (!result.reference_loaded)
    {
        return result;
    }

    const StdVec<Vecd> probe_positions = buildTeam7NativeProbeLinePositionsSi(definition);
    result.samples =
        sampleBProbeFromNativeReloadCase(case_setup, probe_positions, spec.probe_b_sample_air_particles_only);
    if (spec.use_frozen_team7_bz_sign)
    {
        result.sign_audit = buildBzSignConventionAudit(result.samples, result.reference, spec.tesla_to_mT_scale,
                                                       spec.frozen_sign_bz_real, spec.frozen_sign_bz_imag);
        if (spec.run_sign_scan_diagnostic)
        {
            result.sign_scan_audit =
                pickBestBzSignConvention(result.samples, result.reference, spec.tesla_to_mT_scale);
            result.sign_scan_audit_valid = true;
        }
    }
    else
    {
        result.sign_audit = pickBestBzSignConvention(result.samples, result.reference, spec.tesla_to_mT_scale);
    }
    for (const AphiProbeLineSample &sample : result.samples)
    {
        result.max_nearest_distance_mm =
            std::max(result.max_nearest_distance_mm, Real(1.0e3) * sample.neighbor_distance);
    }

    const std::string &report_dir = spec.report_dir;
    std::error_code mkdir_error;
    std::filesystem::create_directories(report_dir, mkdir_error);
    const std::string probe_csv_path = report_dir + "/probe_" + definition.probe_id + "_Bz.csv";
    const std::string comparison_csv_path = report_dir + "/comparison_" + definition.probe_id + "_vs_reference.csv";
    result.probe_csv_written = writeNativeBzProbeCsv(probe_csv_path, result.samples, result.sign_audit);
    result.comparison_csv_written =
        writeNativeBzComparisonCsv(comparison_csv_path, result.samples, result.reference, result.sign_audit);
    return result;
}

inline Real maxAbsSignedBzProbeLineMt(const AphiTeam7NativeBzProbeLineResult &line,
                                      const AphiTeam7BzSignConventionAudit &sign_audit)
{
    Real max_abs = 0.0;
    for (const AphiProbeLineSample &sample : line.samples)
    {
        const Real bz_real_mT = sign_audit.sign_bz_real * sign_audit.tesla_to_mT_scale * sample.b_real[2];
        const Real bz_imag_mT = sign_audit.sign_bz_imag * sign_audit.tesla_to_mT_scale * sample.b_imag[2];
        max_abs = std::max(max_abs, std::max(std::fabs(bz_real_mT), std::fabs(bz_imag_mT)));
    }
    return max_abs;
}

inline AphiTeam7NativeBzReferenceProbeSummary runTeam7NativeBzReferenceProbe(int ac, char *av[],
                                                                             const AphiTeam7NativeBzReferenceProbeSpec &spec)
{
    AphiTeam7NativeBzReferenceProbeSummary summary;
    summary.reload_files_present = allNativeReloadXmlExist();
    summary.reference_data_present =
        std::filesystem::exists("./reference_data/team7/TEAM7_Bz_A1_B1_reference_mT.csv") &&
        std::filesystem::exists("./reference_data/team7/TEAM7_Bz_A2_B2_reference_mT.csv");
    if (!summary.reload_files_present || !summary.reference_data_present)
    {
        return summary;
    }

    AphiTeam7NativeReloadContactCase case_setup(ac, av);
    AphiTeam7NativeBzReferenceProbeSpec effective_spec = spec;
    if (!effective_spec.write_vtp)
    {
        effective_spec.write_vtp = team7NativeCliRequestsWriteVtp(ac, av);
    }
    if (std::abs(effective_spec.max_probe_nearest_distance_m -
                 3.0 * AphiTeam7NativeGeometryConstants::dp_0_default) <= TinyReal)
    {
        effective_spec.max_probe_nearest_distance_m =
            3.0 * team7NativeProbeQuerySpacingM(case_setup.geometry_config);
    }
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    const AphiGMRESResult solver_result =
        runTeam7NativeReloadEmSolve(case_setup, effective_spec.solve_spec, names, joule_names);
    summary.final_true_rel = solver_result.final_true_relative_residual;
    summary.solver_ok = gmresConvergencePassed(solver_result, effective_spec.solve_spec.tolerance) ||
                        summary.final_true_rel <= effective_spec.solve_spec.diagnostic_max_true_rel;

    if (parseTeam7EnvString("TEAM7_SOLVE_AUDIT", "0") == "1")
    {
        AphiLhsAssemblyOptions audit_options;
        audit_options.omega = effective_spec.solve_spec.omega;
        audit_options.use_phi_gauge_penalty = true;
        audit_options.phi_gauge_penalty = effective_spec.solve_spec.phi_gauge_penalty;
        const AphiTeam7NativeAirSolveFieldAudit field_audit =
            auditTeam7NativeAirSolveFields(case_setup, names, audit_options);
        printTeam7NativeAirSolveFieldAudit("post_gmres", field_audit);
    }

    execCurlBOnNativeReloadCase(case_setup, names);

    summary.a1_b1.definition.probe_id = "A1_B1";
    summary.a1_b1.definition.y_m = 0.072;
    summary.a1_b1.definition.z_m = 0.034;
    summary.a1_b1.definition.reference_csv_path = "./reference_data/team7/TEAM7_Bz_A1_B1_reference_mT.csv";
    summary.a2_b2.definition.probe_id = "A2_B2";
    summary.a2_b2.definition.y_m = 0.144;
    summary.a2_b2.definition.z_m = 0.034;
    summary.a2_b2.definition.reference_csv_path = "./reference_data/team7/TEAM7_Bz_A2_B2_reference_mT.csv";

    summary.a1_b1 = runNativeBzProbeLine(case_setup, summary.a1_b1.definition, effective_spec);
    summary.a2_b2 = runNativeBzProbeLine(case_setup, summary.a2_b2.definition, effective_spec);

    if (effective_spec.write_vtp)
    {
        writeTeam7NativeReloadPostSolveVtp(case_setup, names, joule_names);
        summary.vtp_written = true;
    }

    const std::string summary_csv_path = effective_spec.report_dir + "/summary.csv";
    std::ofstream summary_csv(summary_csv_path);
    if (summary_csv)
    {
        summary_csv << std::setprecision(10);
        summary_csv << "metric,value\n";
        summary_csv << "final_true_rel," << summary.final_true_rel << "\n";
        summary_csv << "solver_ok," << (summary.solver_ok ? 1 : 0) << "\n";
        summary_csv << "use_frozen_team7_bz_sign," << (spec.use_frozen_team7_bz_sign ? 1 : 0) << "\n";
        summary_csv << "A1_B1_sign_bz_real," << summary.a1_b1.sign_audit.sign_bz_real << "\n";
        summary_csv << "A1_B1_sign_bz_imag," << summary.a1_b1.sign_audit.sign_bz_imag << "\n";
        summary_csv << "A1_B1_profile_rel_err_real," << summary.a1_b1.sign_audit.profile_rel_err_real << "\n";
        summary_csv << "A1_B1_profile_rel_err_imag," << summary.a1_b1.sign_audit.profile_rel_err_imag << "\n";
        summary_csv << "A1_B1_max_point_rel_err_real," << summary.a1_b1.sign_audit.max_point_rel_err_real << "\n";
        summary_csv << "A1_B1_max_point_rel_err_imag," << summary.a1_b1.sign_audit.max_point_rel_err_imag << "\n";
        summary_csv << "A1_B1_worst_point_rel_err_x_mm," << summary.a1_b1.sign_audit.worst_point_rel_err_x_mm << "\n";
        summary_csv << "A1_B1_max_abs_Bz_mT,"
                    << maxAbsSignedBzProbeLineMt(summary.a1_b1, summary.a1_b1.sign_audit) << "\n";
        summary_csv << "A1_B1_max_nearest_distance_mm," << summary.a1_b1.max_nearest_distance_mm << "\n";
        if (summary.a1_b1.sign_scan_audit_valid)
        {
            summary_csv << "A1_B1_scan_sign_bz_real," << summary.a1_b1.sign_scan_audit.sign_bz_real << "\n";
            summary_csv << "A1_B1_scan_sign_bz_imag," << summary.a1_b1.sign_scan_audit.sign_bz_imag << "\n";
            summary_csv << "A1_B1_scan_profile_rel_err_real," << summary.a1_b1.sign_scan_audit.profile_rel_err_real
                        << "\n";
            summary_csv << "A1_B1_scan_profile_rel_err_imag," << summary.a1_b1.sign_scan_audit.profile_rel_err_imag
                        << "\n";
            summary_csv << "A1_B1_frozen_sign_matches_scan,"
                        << (frozenTeam7BzSignMatchesScan(summary.a1_b1.sign_audit, summary.a1_b1.sign_scan_audit) ? 1
                            : 0)
                        << "\n";
        }
        summary_csv << "A2_B2_profile_rel_err_real," << summary.a2_b2.sign_audit.profile_rel_err_real << "\n";
        summary_csv << "A2_B2_profile_rel_err_imag," << summary.a2_b2.sign_audit.profile_rel_err_imag << "\n";
        summary_csv << "A2_B2_max_nearest_distance_mm," << summary.a2_b2.max_nearest_distance_mm << "\n";
        summary_csv << "vtp_written," << (summary.vtp_written ? 1 : 0) << "\n";
        summary_csv << "vtp_output_dir,output\n";
        summary.summary_csv_written = true;
    }
    return summary;
}

inline bool team7NativeBzReferenceProbePassed(const AphiTeam7NativeBzReferenceProbeSummary &summary,
                                              const AphiTeam7NativeBzReferenceProbeSpec &spec)
{
    const bool a1_samples_ok = summary.a1_b1.reference_loaded && !summary.a1_b1.samples.empty() &&
                               probeSamplesFinite(summary.a1_b1.samples) &&
                               summary.a1_b1.max_nearest_distance_mm <= 1.0e3 * spec.max_probe_nearest_distance_m;
    const bool a1_profile_ok =
        !spec.require_a1b1_profile_gate ||
        (summary.a1_b1.sign_audit.profile_rel_err_real <= spec.max_a1b1_profile_rel_err &&
         summary.a1_b1.sign_audit.profile_rel_err_imag <= spec.max_a1b1_profile_rel_err);
    const Real max_abs_bz_a1b1_mT = maxAbsSignedBzProbeLineMt(summary.a1_b1, summary.a1_b1.sign_audit);
    const bool bz_magnitude_ok =
        !spec.require_bz_magnitude_sanity ||
        (max_abs_bz_a1b1_mT > spec.min_max_abs_bz_a1b1_mT && max_abs_bz_a1b1_mT < spec.max_max_abs_bz_a1b1_mT);
    const bool sign_frozen_ok =
        !spec.use_frozen_team7_bz_sign || !spec.require_frozen_sign_matches_scan ||
        !summary.a1_b1.sign_scan_audit_valid ||
        frozenTeam7BzSignMatchesScan(summary.a1_b1.sign_audit, summary.a1_b1.sign_scan_audit);
    const bool csv_ok = summary.a1_b1.probe_csv_written && summary.a1_b1.comparison_csv_written &&
                        summary.a2_b2.probe_csv_written && summary.a2_b2.comparison_csv_written &&
                        summary.summary_csv_written;
    return summary.reload_files_present && summary.reference_data_present && summary.solver_ok && a1_samples_ok &&
           a1_profile_ok && bz_magnitude_ok && sign_frozen_ok && csv_ok;
}

inline void printTeam7NativeBzReferenceProbeSummary(const char *test_name,
                                                    const AphiTeam7NativeBzReferenceProbeSummary &summary,
                                                    const AphiTeam7NativeBzReferenceProbeSpec &spec, bool passed)
{
    std::cout << test_name << " passed=" << (passed ? 1 : 0)
              << " reload_files_present=" << (summary.reload_files_present ? 1 : 0)
              << " reference_data_present=" << (summary.reference_data_present ? 1 : 0)
              << " solver_ok=" << (summary.solver_ok ? 1 : 0) << " final_true_rel=" << summary.final_true_rel
              << " frozen_sign=" << (spec.use_frozen_team7_bz_sign ? 1 : 0)
              << " A1_B1_sign=(" << summary.a1_b1.sign_audit.sign_bz_real << "," << summary.a1_b1.sign_audit.sign_bz_imag
              << ") A1_B1_profile_rel_err_real=" << summary.a1_b1.sign_audit.profile_rel_err_real
              << " A1_B1_profile_rel_err_imag=" << summary.a1_b1.sign_audit.profile_rel_err_imag
              << " A1_B1_max_point_rel_err_real=" << summary.a1_b1.sign_audit.max_point_rel_err_real
              << " worst_x_mm=" << summary.a1_b1.sign_audit.worst_point_rel_err_x_mm
              << " A1_B1_max_nearest_distance_mm=" << summary.a1_b1.max_nearest_distance_mm
              << " A2_B2_profile_rel_err_real=" << summary.a2_b2.sign_audit.profile_rel_err_real
              << " A2_B2_profile_rel_err_imag=" << summary.a2_b2.sign_audit.profile_rel_err_imag
              << " A1_B1_max_abs_Bz_mT=" << maxAbsSignedBzProbeLineMt(summary.a1_b1, summary.a1_b1.sign_audit);
    if (summary.a1_b1.sign_scan_audit_valid)
    {
        std::cout << " A1_B1_scan_sign=(" << summary.a1_b1.sign_scan_audit.sign_bz_real << ","
                  << summary.a1_b1.sign_scan_audit.sign_bz_imag
                  << ") frozen_matches_scan="
                  << (frozenTeam7BzSignMatchesScan(summary.a1_b1.sign_audit, summary.a1_b1.sign_scan_audit) ? 1 : 0);
    }
    std::cout << " vtp_written=" << (summary.vtp_written ? 1 : 0);
    if (summary.vtp_written)
    {
        std::cout << " vtp_output_dir=output";
    }
    std::cout << " report_dir=" << spec.report_dir << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_NATIVE_BZ_REFERENCE_PROBE_HELPERS_H
