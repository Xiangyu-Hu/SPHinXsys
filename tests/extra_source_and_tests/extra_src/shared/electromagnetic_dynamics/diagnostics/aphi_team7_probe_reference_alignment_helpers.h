#ifndef APHI_TEAM7_PROBE_REFERENCE_ALIGNMENT_HELPERS_H
#define APHI_TEAM7_PROBE_REFERENCE_ALIGNMENT_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_team7_canonical_case_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiTeam7ProbeReferenceAlignmentSummary
{
    std::string report_dir = "team7_probe_reference_alignment_report";
    Real reference_dp = 0.0;
    Real candidate_dp = 0.0;
    Real centerline_yz_band_reference = 0.0;
    Real centerline_yz_band_candidate = 0.0;
    /** Canonical quantitative path: binned conductor A-phi block norm along x. */
    Real centerline_profile_rel_A_block = 0.0;
    /** Probe nearest-particle path: B_abs along physical-box center (y,z at box mid). */
    Real probe_box_center_profile_rel_B_abs = 0.0;
    Real probe_box_center_profile_rel_Bz_real = 0.0;
    Real probe_box_center_profile_rel_Bz_imag = 0.0;
    /** Probe path aligned with canonical conductor mid-plane (y_mid,z_mid of conductor box). */
    Real probe_conductor_midplane_profile_rel_B_abs = 0.0;
    Real probe_conductor_midplane_profile_rel_Bz_real = 0.0;
    bool reference_converged = false;
    bool candidate_converged = false;
    bool metadata_csv_written = false;
    bool aligned_csv_written = false;
    bool box_center_probe_csv_written = false;
    bool summary_csv_written = false;
    bool all_values_finite = false;
    bool diagnostic_only = true;
};

inline AphiSourceDrivenEmSolveSpec team7CanonicalSourceDrivenSpec(Real dp)
{
    using benchmark::AphiTeam7CanonicalCaseSpec;
    AphiSourceDrivenEmSolveSpec spec;
    spec.dp = dp;
    spec.omega = AphiTeam7CanonicalCaseSpec::omega;
    spec.phi_gauge_penalty = AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    spec.impressed_current_amplitude = AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    spec.tolerance = AphiTeam7CanonicalCaseSpec::tolerance;
    spec.restart_dimension = AphiTeam7CanonicalCaseSpec::restart_dimension;
    spec.max_outer_iterations = AphiTeam7CanonicalCaseSpec::max_outer_iterations;
    spec.write_vtp = false;
    spec.write_probe_csv = false;
    spec.max_air_to_conductor_joule_ratio = 1.0;
    return spec;
}

inline StdVec<Vecd> buildTeam7ConductorMidplaneCenterlinePositions(
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, UnsignedInt sample_count)
{
    const Real y_mid = Real(0.5) * (layout.conductor.ymin + layout.conductor.ymax);
    const Real z_mid = Real(0.5) * (layout.conductor.zmin + layout.conductor.zmax);
    const Real x_min = layout.conductor.xmin;
    const Real x_max = layout.conductor.xmax;
    StdVec<Vecd> positions;
    positions.reserve(static_cast<size_t>(sample_count));
    for (UnsignedInt sample_index = 0; sample_index < sample_count; ++sample_index)
    {
        const Real alpha =
            static_cast<Real>(sample_index) / (static_cast<Real>(sample_count - 1) + TinyReal);
        positions.push_back(Vecd(x_min + alpha * (x_max - x_min), y_mid, z_mid));
    }
    return positions;
}

inline Real probeLineProfileL2RelativeDifference(
    const StdVec<AphiProbeLineSample> &reference, const StdVec<AphiProbeLineSample> &approximate,
    const std::function<Real(const AphiProbeLineSample &)> &metric)
{
    if (reference.size() != approximate.size() || reference.empty())
    {
        return Real(0);
    }
    Real sum_sq = 0.0;
    Real ref_sum_sq = 0.0;
    for (size_t i = 0; i != reference.size(); ++i)
    {
        const Real ref_value = metric(reference[i]);
        const Real approx_value = metric(approximate[i]);
        if (!std::isfinite(ref_value) || !std::isfinite(approx_value))
        {
            return std::numeric_limits<Real>::infinity();
        }
        const Real diff = approx_value - ref_value;
        sum_sq += diff * diff;
        ref_sum_sq += ref_value * ref_value;
    }
    return std::sqrt(sum_sq) / (std::sqrt(ref_sum_sq) + TinyReal);
}

inline StdVec<AphiProbeLineSample> runTeam7SourceDrivenProbeLineSamples(
    int ac, char *av[], const AphiSourceDrivenEmSolveSpec &spec,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, const StdVec<Vecd> &sample_positions)
{
    const Real boundary_width = 3.0 * spec.dp * spec.boundary_width_scale;
    AphiJouleHeatingFieldNames joule_fields;
    AphiSourceDrivenEmSolveFieldNames obs_fields;
    AphiVariableNames names;

    AphiLhsTestBody test_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    (void)register_joule_fields;
    (void)execSourceDrivenEmJoulePipelineOnBody(test_body, layout, spec, names, joule_fields, &obs_fields);

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    AphiProbeObservedFieldNames probe_fields;
    probe_fields.b_real = obs_fields.b_real;
    probe_fields.b_imag = obs_fields.b_imag;
    probe_fields.b_magnitude = obs_fields.b_magnitude;
    probe_fields.material_region_id = obs_fields.material_region_id;
    return sampleLineProbesNearestParticle(particles, positions, total_real_particles, sample_positions, probe_fields,
                                           joule_fields, names.material);
}

inline bool writeTeam7ProbeReferenceMetadataCsv(const std::string &path,
                                                const AphiTeam7ProbeReferenceAlignmentSummary &summary)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "key,value\n";
    output << "reference_dp," << summary.reference_dp << "\n";
    output << "candidate_dp," << summary.candidate_dp << "\n";
    output << "centerline_yz_band_reference," << summary.centerline_yz_band_reference << "\n";
    output << "centerline_yz_band_candidate," << summary.centerline_yz_band_candidate << "\n";
    output << "canonical_metric,binned_conductor_Aphi_block_norm\n";
    output << "probe_box_center_metric,nearest_particle_B_from_curlA\n";
    output << "probe_conductor_midplane_metric,nearest_particle_B_at_conductor_yz_mid\n";
    output << "aligned_csv_x_axis,conductor_x_bin_centers\n";
    output << "box_center_probe_csv,aphi_team7_probe_reference_box_center_probe.csv\n";
    output << "probe_sampling,nearest_real_particle\n";
    output << "b_postprocess,BCorrectedGrad_curlA\n";
    output << "frequency_omega," << benchmark::AphiTeam7CanonicalCaseSpec::omega << "\n";
    output << "impressed_current_amplitude," << benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude
           << "\n";
    output << "diagnostic_only," << (summary.diagnostic_only ? 1 : 0) << "\n";
    output << "standard_team7_validation,0\n";
    return true;
}

inline bool writeTeam7ProbeReferenceAlignedCsv(const std::string &path,
                                              const AphiTeam7CenterlineProfile &reference_centerline,
                                              const AphiTeam7CenterlineProfile &candidate_centerline,
                                              const StdVec<AphiProbeLineSample> &reference_midplane_probe,
                                              const StdVec<AphiProbeLineSample> &candidate_midplane_probe)
{
    const size_t row_count = reference_centerline.x_centers.size();
    if (row_count == 0 || candidate_centerline.x_centers.size() != row_count ||
        reference_midplane_probe.size() != candidate_midplane_probe.size() ||
        reference_midplane_probe.size() != row_count)
    {
        return false;
    }

    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "x_conductor_midplane,reference_A_block_norm,candidate_A_block_norm,reference_probe_midplane_B_abs,"
              "candidate_probe_midplane_B_abs,reference_probe_midplane_Bz_real,candidate_probe_midplane_Bz_real,"
              "reference_probe_midplane_Bz_imag,candidate_probe_midplane_Bz_imag,reference_probe_nearest_distance,"
              "candidate_probe_nearest_distance\n";
    for (size_t i = 0; i != row_count; ++i)
    {
        output << reference_centerline.x_centers[i] << "," << reference_centerline.field_norm[i] << ","
               << candidate_centerline.field_norm[i] << "," << reference_midplane_probe[i].b_abs << ","
               << candidate_midplane_probe[i].b_abs << "," << reference_midplane_probe[i].b_real[2] << ","
               << candidate_midplane_probe[i].b_real[2] << "," << reference_midplane_probe[i].b_imag[2] << ","
               << candidate_midplane_probe[i].b_imag[2] << "," << reference_midplane_probe[i].neighbor_distance << ","
               << candidate_midplane_probe[i].neighbor_distance << "\n";
    }
    return true;
}

inline bool writeTeam7ProbeReferenceBoxCenterProbeCsv(const std::string &path,
                                                      const StdVec<AphiProbeLineSample> &reference_box_probe,
                                                      const StdVec<AphiProbeLineSample> &candidate_box_probe)
{
    if (reference_box_probe.size() != candidate_box_probe.size() || reference_box_probe.empty())
    {
        return false;
    }
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "x_box_center,y_box_center,z_box_center,reference_B_abs,candidate_B_abs,reference_Bz_real,"
              "candidate_Bz_real,reference_nearest_distance,candidate_nearest_distance\n";
    for (size_t i = 0; i != reference_box_probe.size(); ++i)
    {
        output << reference_box_probe[i].position[0] << "," << reference_box_probe[i].position[1] << ","
               << reference_box_probe[i].position[2] << "," << reference_box_probe[i].b_abs << ","
               << candidate_box_probe[i].b_abs << "," << reference_box_probe[i].b_real[2] << ","
               << candidate_box_probe[i].b_real[2] << "," << reference_box_probe[i].neighbor_distance << ","
               << candidate_box_probe[i].neighbor_distance << "\n";
    }
    return true;
}

inline bool writeTeam7ProbeReferenceProfileRelSummaryCsv(const std::string &path,
                                                         const AphiTeam7ProbeReferenceAlignmentSummary &summary)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "metric,value\n";
    output << "centerline_profile_rel_A_block," << summary.centerline_profile_rel_A_block << "\n";
    output << "probe_box_center_profile_rel_B_abs," << summary.probe_box_center_profile_rel_B_abs << "\n";
    output << "probe_box_center_profile_rel_Bz_real," << summary.probe_box_center_profile_rel_Bz_real << "\n";
    output << "probe_box_center_profile_rel_Bz_imag," << summary.probe_box_center_profile_rel_Bz_imag << "\n";
    output << "probe_conductor_midplane_profile_rel_B_abs," << summary.probe_conductor_midplane_profile_rel_B_abs
           << "\n";
    output << "probe_conductor_midplane_profile_rel_Bz_real," << summary.probe_conductor_midplane_profile_rel_Bz_real
           << "\n";
    output << "reference_converged," << (summary.reference_converged ? 1 : 0) << "\n";
    output << "candidate_converged," << (summary.candidate_converged ? 1 : 0) << "\n";
    return true;
}

inline bool team7ProbeReferenceAlignmentFinite(const AphiTeam7ProbeReferenceAlignmentSummary &summary)
{
    const auto finite = [](Real value) { return std::isfinite(value); };
    return finite(summary.centerline_profile_rel_A_block) && finite(summary.probe_box_center_profile_rel_B_abs) &&
           finite(summary.probe_box_center_profile_rel_Bz_real) &&
           finite(summary.probe_box_center_profile_rel_Bz_imag) &&
           finite(summary.probe_conductor_midplane_profile_rel_B_abs) &&
           finite(summary.probe_conductor_midplane_profile_rel_Bz_real);
}

inline AphiTeam7ProbeReferenceAlignmentSummary runTeam7ProbeReferenceAlignmentAudit(
    int ac, char *av[], const std::string &report_dir = "team7_probe_reference_alignment_report")
{
    using benchmark::AphiTeam7CanonicalCaseSpec;
    std::error_code mkdir_error;
    std::filesystem::create_directories(report_dir, mkdir_error);

    AphiTeam7ProbeReferenceAlignmentSummary summary;
    summary.report_dir = report_dir;
    summary.reference_dp = AphiTeam7CanonicalCaseSpec::reference_dp;
    summary.candidate_dp = AphiTeam7CanonicalCaseSpec::candidate_dp;
    summary.centerline_yz_band_reference = AphiTeam7CanonicalCaseSpec::centerline_yz_band_factor * summary.reference_dp;
    summary.centerline_yz_band_candidate = AphiTeam7CanonicalCaseSpec::centerline_yz_band_factor * summary.candidate_dp;

    const AphiTeam7CanonicalCaseRunResult reference_run =
        runTeam7PhysicalCanonicalCase(ac, av, summary.reference_dp);
    const AphiTeam7CanonicalCaseRunResult candidate_run = runTeam7PhysicalCanonicalCase(ac, av, summary.candidate_dp);
    summary.reference_converged = reference_run.metrics.converged;
    summary.candidate_converged = candidate_run.metrics.converged;
    summary.centerline_profile_rel_A_block =
        hostCenterlineProfileL2RelativeDifference(reference_run.centerline, candidate_run.centerline);

    const AphiSourceDrivenEmSolveSpec reference_spec = team7CanonicalSourceDrivenSpec(summary.reference_dp);
    const AphiSourceDrivenEmSolveSpec candidate_spec = team7CanonicalSourceDrivenSpec(summary.candidate_dp);
    const benchmark::AphiTeam7LikeUnitBoxLayout layout = reference_run.layout;

    const UnsignedInt probe_sample_count = AphiTeam7CanonicalCaseSpec::centerline_bins;
    const StdVec<Vecd> box_center_positions = buildTeam7CenterlineSamplePositions(
        AphiTeam7CanonicalCaseSpec::body_length, AphiTeam7CanonicalCaseSpec::body_height,
        AphiTeam7CanonicalCaseSpec::body_width, probe_sample_count);
    const StdVec<Vecd> conductor_midplane_positions =
        buildTeam7ConductorMidplaneCenterlinePositions(layout, probe_sample_count);

    const StdVec<AphiProbeLineSample> reference_box_probe =
        runTeam7SourceDrivenProbeLineSamples(ac, av, reference_spec, layout, box_center_positions);
    const StdVec<AphiProbeLineSample> candidate_box_probe =
        runTeam7SourceDrivenProbeLineSamples(ac, av, candidate_spec, layout, box_center_positions);
    const StdVec<AphiProbeLineSample> reference_midplane_probe =
        runTeam7SourceDrivenProbeLineSamples(ac, av, reference_spec, layout, conductor_midplane_positions);
    const StdVec<AphiProbeLineSample> candidate_midplane_probe =
        runTeam7SourceDrivenProbeLineSamples(ac, av, candidate_spec, layout, conductor_midplane_positions);

    summary.probe_box_center_profile_rel_B_abs =
        probeLineProfileL2RelativeDifference(reference_box_probe, candidate_box_probe,
                                             [](const AphiProbeLineSample &sample) { return sample.b_abs; });
    summary.probe_box_center_profile_rel_Bz_real =
        probeLineProfileL2RelativeDifference(reference_box_probe, candidate_box_probe,
                                             [](const AphiProbeLineSample &sample) { return sample.b_real[2]; });
    summary.probe_box_center_profile_rel_Bz_imag =
        probeLineProfileL2RelativeDifference(reference_box_probe, candidate_box_probe,
                                             [](const AphiProbeLineSample &sample) { return sample.b_imag[2]; });
    summary.probe_conductor_midplane_profile_rel_B_abs =
        probeLineProfileL2RelativeDifference(reference_midplane_probe, candidate_midplane_probe,
                                             [](const AphiProbeLineSample &sample) { return sample.b_abs; });
    summary.probe_conductor_midplane_profile_rel_Bz_real =
        probeLineProfileL2RelativeDifference(reference_midplane_probe, candidate_midplane_probe,
                                             [](const AphiProbeLineSample &sample) { return sample.b_real[2]; });

    summary.metadata_csv_written =
        writeTeam7ProbeReferenceMetadataCsv(report_dir + "/aphi_team7_probe_reference_metadata.csv", summary);
    summary.aligned_csv_written = writeTeam7ProbeReferenceAlignedCsv(
        report_dir + "/aphi_team7_probe_reference_aligned_comparison.csv", reference_run.centerline,
        candidate_run.centerline, reference_midplane_probe, candidate_midplane_probe);
    summary.box_center_probe_csv_written = writeTeam7ProbeReferenceBoxCenterProbeCsv(
        report_dir + "/aphi_team7_probe_reference_box_center_probe.csv", reference_box_probe, candidate_box_probe);
    summary.summary_csv_written = writeTeam7ProbeReferenceProfileRelSummaryCsv(
        report_dir + "/aphi_team7_probe_reference_profile_rel_summary.csv", summary);
    summary.all_values_finite = team7ProbeReferenceAlignmentFinite(summary);
    return summary;
}

inline bool team7ProbeReferenceAlignmentPassed(const AphiTeam7ProbeReferenceAlignmentSummary &summary)
{
    return summary.reference_converged && summary.candidate_converged && summary.metadata_csv_written &&
           summary.aligned_csv_written && summary.box_center_probe_csv_written && summary.summary_csv_written &&
           summary.all_values_finite;
}

inline void printTeam7ProbeReferenceAlignmentSummary(const char *test_name,
                                                     const AphiTeam7ProbeReferenceAlignmentSummary &summary)
{
    std::cout << test_name << " reference_dp=" << summary.reference_dp << " candidate_dp=" << summary.candidate_dp
              << " centerline_profile_rel_A_block=" << summary.centerline_profile_rel_A_block
              << " probe_box_center_profile_rel_B_abs=" << summary.probe_box_center_profile_rel_B_abs
              << " probe_box_center_profile_rel_Bz_real=" << summary.probe_box_center_profile_rel_Bz_real
              << " probe_conductor_midplane_profile_rel_B_abs=" << summary.probe_conductor_midplane_profile_rel_B_abs
              << " probe_conductor_midplane_profile_rel_Bz_real=" << summary.probe_conductor_midplane_profile_rel_Bz_real
              << " metadata_csv_written=" << (summary.metadata_csv_written ? 1 : 0)
              << " aligned_csv_written=" << (summary.aligned_csv_written ? 1 : 0)
              << " box_center_probe_csv_written=" << (summary.box_center_probe_csv_written ? 1 : 0)
              << " summary_csv_written=" << (summary.summary_csv_written ? 1 : 0)
              << " all_values_finite=" << (summary.all_values_finite ? 1 : 0)
              << " diagnostic_only=" << (summary.diagnostic_only ? 1 : 0) << " report_dir=" << summary.report_dir
              << " standard_team7_validation=0" << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
