#ifndef APHI_TEAM7_SINGLE_BODY_VS_CONTACT_COMPARISON_HELPERS_H
#define APHI_TEAM7_SINGLE_BODY_VS_CONTACT_COMPARISON_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_team7_canonical_case_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_probe_metric_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_probe_reference_alignment_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_physical_region_audit_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiTeam7SingleBodyVsContactComparisonSpec
{
    Real dp = 0.1;
    Real body_length = benchmark::AphiTeam7PhysicalDimensions::length;
    Real body_height = benchmark::AphiTeam7PhysicalDimensions::height;
    Real body_width = benchmark::AphiTeam7PhysicalDimensions::width;
    Real max_conductor_joule_rel = 0.20;
    Real max_b_probe_profile_rel = 0.30;
    UnsignedInt probe_sample_count = benchmark::AphiTeam7CanonicalCaseSpec::centerline_bins;
    std::string report_dir = "team7_single_body_vs_contact_report";
};

struct AphiTeam7SingleBodyVsContactRunSnapshot
{
    bool converged = false;
    UnsignedInt num_iterations = 0;
    Real final_residual = 0.0;
    Real conductor_joule_integral = 0.0;
    Real air_joule_integral = 0.0;
    Real source_rhs_l2 = 0.0;
    Real max_abs_J = 0.0;
    bool finite_fields = false;
    AphiConductorInterfaceSpikeMetrics interface_spike{};
    StdVec<AphiProbeLineSample> conductor_midplane_b_probe{};
};

struct AphiTeam7SingleBodyVsContactComparisonSummary
{
    AphiTeam7SingleBodyVsContactComparisonSpec spec{};
    AphiTeam7SingleBodyVsContactRunSnapshot single_body{};
    AphiTeam7SingleBodyVsContactRunSnapshot three_body_contact{};
    AphiConductorInterfaceSpikeMetrics interface_spike{};
    Real conductor_joule_rel = 0.0;
    Real b_probe_profile_rel = 0.0;
    bool comparison_csv_written = false;
    bool summary_csv_written = false;
};

inline AphiSourceDrivenEmSolveSpec team7ComparisonSourceDrivenSpec(const AphiTeam7SingleBodyVsContactComparisonSpec &spec)
{
    AphiSourceDrivenEmSolveSpec solve_spec;
    solve_spec.dp = spec.dp;
    solve_spec.body_length = spec.body_length;
    solve_spec.body_height = spec.body_height;
    solve_spec.body_width = spec.body_width;
    solve_spec.omega = benchmark::AphiTeam7CanonicalCaseSpec::omega;
    solve_spec.phi_gauge_penalty = benchmark::AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    solve_spec.impressed_current_amplitude = benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    solve_spec.tolerance = benchmark::AphiTeam7CanonicalCaseSpec::tolerance;
    solve_spec.restart_dimension = benchmark::AphiTeam7CanonicalCaseSpec::restart_dimension;
    solve_spec.max_outer_iterations = benchmark::AphiTeam7CanonicalCaseSpec::max_outer_iterations;
    solve_spec.write_vtp = false;
    solve_spec.write_probe_csv = false;
    solve_spec.max_air_to_conductor_joule_ratio = 1.0;
    return solve_spec;
}

inline StdVec<AphiProbeLineSample> sampleBProbeLineNearestParticle(
    BaseParticles &particles, const Vecd *positions, size_t total_real_particles,
    const StdVec<Vecd> &sample_positions, const std::string &b_real_name, const std::string &b_imag_name,
    const std::function<bool(const Vecd &)> &particle_filter = [](const Vecd &) { return true; })
{
    syncVariableToHost<Vecd>(particles, b_real_name);
    syncVariableToHost<Vecd>(particles, b_imag_name);
    const Vecd *b_real = particles.getVariableDataByName<Vecd>(b_real_name);
    const Vecd *b_imag = particles.getVariableDataByName<Vecd>(b_imag_name);

    StdVec<AphiProbeLineSample> samples;
    samples.reserve(sample_positions.size());
    for (const Vecd &sample_position : sample_positions)
    {
        size_t index = 0;
        Real nearest_distance_squared = std::numeric_limits<Real>::max();
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            if (!particle_filter(positions[i]))
            {
                continue;
            }
            const Real distance_squared = (positions[i] - sample_position).squaredNorm();
            if (distance_squared < nearest_distance_squared)
            {
                nearest_distance_squared = distance_squared;
                index = i;
            }
        }
        AphiProbeLineSample sample;
        sample.position = sample_position;
        sample.nearest_particle_index = index;
        sample.neighbor_distance = (positions[index] - sample_position).norm();
        sample.b_real = b_real[index];
        sample.b_imag = b_imag[index];
        sample.b_abs = vec3Abs(sample.b_real);
        samples.push_back(sample);
    }
    return samples;
}

inline void execCurlBOnThreeBodyCase(AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names)
{
    static constexpr const char *k_b_real = "ComparisonBReal";
    static constexpr const char *k_b_imag = "ComparisonBImag";
    execBodyCurlBFromADiagnostic(case_setup.air_body(), case_setup.air_inner(), names, k_b_real, k_b_imag,
                                 AphiBCurlDiagnosticMode::BCorrectedGrad);
    execBodyCurlBFromADiagnostic(case_setup.coil_body, case_setup.coil_inner(), names, k_b_real, k_b_imag,
                                 AphiBCurlDiagnosticMode::BCorrectedGrad);
    execBodyCurlBFromADiagnostic(case_setup.plate_body, case_setup.plate_inner(), names, k_b_real, k_b_imag,
                                 AphiBCurlDiagnosticMode::BCorrectedGrad);
}

inline StdVec<AphiProbeLineSample> sampleBProbeFromThreeBodyCase(AphiTeam7ThreeBodyContactCase &case_setup,
                                                                 const StdVec<Vecd> &sample_positions)
{
    static constexpr const char *k_b_real = "ComparisonBReal";
    static constexpr const char *k_b_imag = "ComparisonBImag";
    BaseParticles *const views[] = {&case_setup.air_body().getBaseParticles(), &case_setup.coil_body.getBaseParticles(),
                                    &case_setup.plate_body.getBaseParticles()};

    StdVec<AphiProbeLineSample> samples;
    samples.reserve(sample_positions.size());
    for (const Vecd &sample_position : sample_positions)
    {
        AphiProbeLineSample sample;
        sample.position = sample_position;
        Real best_distance_squared = std::numeric_limits<Real>::max();
        for (BaseParticles *particles : views)
        {
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

inline AphiTeam7SingleBodyVsContactRunSnapshot runSingleBodyTeam7ComparisonSnapshot(
    int ac, char *av[], const AphiTeam7SingleBodyVsContactComparisonSpec &spec,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, const StdVec<Vecd> &probe_positions)
{
    const AphiSourceDrivenEmSolveSpec solve_spec = team7ComparisonSourceDrivenSpec(spec);
    const Real boundary_width = 3.0 * spec.dp;
    AphiJouleHeatingFieldNames joule_fields;
    AphiVariableNames names;
    AphiSourceDrivenEmSolveFieldNames obs_fields;
    AphiLhsTestBody test_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    (void)register_joule_fields;
    const AphiMatrixFreeSolverResult em_result =
        execSourceDrivenEmJoulePipelineOnBody(test_body, layout, solve_spec, names, joule_fields, &obs_fields);
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();
    const auto in_physical_conductor = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Conductor);
    };
    const auto in_physical_air = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Air);
    };
    const auto in_physical_source = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Coil);
    };

    AphiTeam7SingleBodyVsContactRunSnapshot snapshot;
    snapshot.converged = gmresConvergencePassed(em_result, solve_spec.tolerance);
    snapshot.num_iterations = em_result.outer_iteration_count;
    snapshot.final_residual = em_result.final_true_relative_residual;
    snapshot.conductor_joule_integral = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_conductor, joule_fields.joule_heat_source);
    snapshot.air_joule_integral = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles,
                                                                        in_physical_air, joule_fields.joule_heat_source);
    snapshot.source_rhs_l2 =
        hostParticleRegionBlockNorm(particles, names.rhs, positions, total_real_particles, in_physical_source);
    snapshot.max_abs_J = hostParticleRegionScalarMax(particles, obs_fields.j_magnitude, positions, total_real_particles,
                                                     in_physical_conductor);
    snapshot.finite_fields =
        hostSourceDrivenFieldsFinite(particles, names, joule_fields, obs_fields, total_real_particles);
    const auto in_conductor = [&](const Vecd &position) {
        return team7ParticleInRegion(position, layout, AphiBenchmarkMaterialRegion::Conductor);
    };
    snapshot.conductor_midplane_b_probe =
        sampleBProbeLineNearestParticle(particles, positions, total_real_particles, probe_positions,
                                        obs_fields.b_real, obs_fields.b_imag, in_conductor);
    return snapshot;
}

inline AphiTeam7SingleBodyVsContactRunSnapshot runThreeBodyContactTeam7ComparisonSnapshot(
    int ac, char *av[], const AphiTeam7SingleBodyVsContactComparisonSpec &spec,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, const StdVec<Vecd> &probe_positions)
{
    const Real boundary_width = 3.0 * spec.dp;
    const Real core_shell = 2.0 * spec.dp;
    const Real omega = benchmark::AphiTeam7CanonicalCaseSpec::omega;
    const Real phi_gauge_penalty = benchmark::AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    const Real tolerance = benchmark::AphiTeam7CanonicalCaseSpec::tolerance;
    const Real impressed_current_amplitude = benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    AphiTeam7ThreeBodyContactCase case_setup(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width,
                                             ac, av);
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = false;

    const AphiMatrixFreeSolverOptions solver_options =
        defaultCoupledContactGMRESOptions(tolerance, benchmark::AphiTeam7CanonicalCaseSpec::restart_dimension,
                                          benchmark::AphiTeam7CanonicalCaseSpec::max_outer_iterations);
    const AphiTeam7ContactCoupledGmresResult result = runTeam7ThreeBodyCoupledContactGmres(
        case_setup, names, joule_names, options, solver_options, spec.body_length, spec.body_height, spec.body_width,
        core_shell, coil_current_real, coil_current_imag, impressed_current_amplitude, true);

    const AphiThreeBodyContactAuditSummary audit =
        hostThreeBodyContactAuditSummary(case_setup, names, joule_names, spec.body_length, spec.body_height, spec.body_width);

    AphiTeam7SingleBodyVsContactRunSnapshot snapshot;
    snapshot.converged = gmresConvergencePassed(result.solver_result, tolerance);
    snapshot.num_iterations = result.solver_result.outer_iteration_count;
    snapshot.final_residual = result.solver_result.final_true_relative_residual;
    snapshot.conductor_joule_integral = result.plate.plate_joule_power;
    snapshot.air_joule_integral = audit.global_air_joule;
    snapshot.source_rhs_l2 = audit.coil.rhs_l2;
    snapshot.max_abs_J = result.plate.plate_j_L2;
    snapshot.finite_fields = std::isfinite(snapshot.conductor_joule_integral) && snapshot.conductor_joule_integral > 0.0;
    execCurlBOnThreeBodyCase(case_setup, names);
    snapshot.conductor_midplane_b_probe = sampleBProbeFromThreeBodyCase(case_setup, probe_positions);
    syncVariableToHost<Vecd>(case_setup.plate_body.getBaseParticles(), "Position");
    const Vecd *plate_positions = case_setup.plate_body.getBaseParticles().getVariableDataByName<Vecd>("Position");
    const size_t plate_count = case_setup.plate_body.getBaseParticles().TotalRealParticles();
    snapshot.interface_spike = hostConductorInterfaceSpikeMetrics(
        case_setup.plate_body.getBaseParticles(), joule_names, plate_positions, plate_count, layout, spec.body_length,
        spec.body_height, spec.body_width, core_shell, layout.conductor.xmin, spec.dp, spec.dp);
    return snapshot;
}

inline bool writeTeam7SingleBodyVsContactComparisonCsv(const std::string &path,
                                                       const AphiTeam7SingleBodyVsContactComparisonSummary &summary)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "case,converged,num_iterations,final_residual,conductor_joule,air_joule,source_rhs_l2,max_abs_J,finite_"
              "fields\n";
    const auto write_row = [&](const char *label, const AphiTeam7SingleBodyVsContactRunSnapshot &run) {
        output << label << "," << (run.converged ? 1 : 0) << "," << run.num_iterations << "," << run.final_residual
               << "," << run.conductor_joule_integral << "," << run.air_joule_integral << "," << run.source_rhs_l2
               << "," << run.max_abs_J << "," << (run.finite_fields ? 1 : 0) << "\n";
    };
    write_row("single_body", summary.single_body);
    write_row("three_body_contact", summary.three_body_contact);
    return true;
}

inline bool writeTeam7SingleBodyVsContactSummaryCsv(const std::string &path,
                                                    const AphiTeam7SingleBodyVsContactComparisonSummary &summary)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "metric,value\n";
    output << "conductor_joule_rel," << summary.conductor_joule_rel << "\n";
    output << "b_probe_profile_rel," << summary.b_probe_profile_rel << "\n";
    output << "max_conductor_joule_rel_gate," << summary.spec.max_conductor_joule_rel << "\n";
    output << "max_b_probe_profile_rel_gate," << summary.spec.max_b_probe_profile_rel << "\n";
    output << "interface_spike_hard_pass," << (summary.interface_spike.spike_hard_pass ? 1 : 0) << "\n";
    output << "lambda_A," << 0 << "\n";
    return true;
}

inline AphiTeam7SingleBodyVsContactComparisonSummary runTeam7SingleBodyVsThreeBodyContactComparison(
    int ac, char *av[], const AphiTeam7SingleBodyVsContactComparisonSpec &spec = {})
{
    std::error_code mkdir_error;
    std::filesystem::create_directories(spec.report_dir, mkdir_error);

    const benchmark::AphiTeam7LikeUnitBoxLayout layout =
        benchmark::buildTeam7LayoutForBox(spec.body_length, spec.body_height, spec.body_width);
    const StdVec<Vecd> probe_positions =
        buildTeam7ConductorMidplaneCenterlinePositions(layout, spec.probe_sample_count);

    AphiTeam7SingleBodyVsContactComparisonSummary summary;
    summary.spec = spec;
    summary.single_body = runSingleBodyTeam7ComparisonSnapshot(ac, av, spec, layout, probe_positions);
    summary.three_body_contact = runThreeBodyContactTeam7ComparisonSnapshot(ac, av, spec, layout, probe_positions);

    summary.conductor_joule_rel =
        relativeMetricChange(summary.single_body.conductor_joule_integral,
                             summary.three_body_contact.conductor_joule_integral);
    summary.b_probe_profile_rel = probeLineProfileL2RelativeDifference(
        summary.single_body.conductor_midplane_b_probe, summary.three_body_contact.conductor_midplane_b_probe,
        [](const AphiProbeLineSample &sample) { return sample.b_abs; });
    summary.interface_spike = summary.three_body_contact.interface_spike;

    summary.comparison_csv_written = writeTeam7SingleBodyVsContactComparisonCsv(
        spec.report_dir + "/aphi_team7_single_body_vs_contact_comparison.csv", summary);
    summary.summary_csv_written =
        writeTeam7SingleBodyVsContactSummaryCsv(spec.report_dir + "/aphi_team7_single_body_vs_contact_summary.csv", summary);
    return summary;
}

inline bool team7SingleBodyVsThreeBodyContactComparisonPassed(
    const AphiTeam7SingleBodyVsContactComparisonSummary &summary)
{
    return summary.single_body.converged && summary.three_body_contact.converged && summary.single_body.finite_fields &&
           summary.three_body_contact.finite_fields && summary.single_body.conductor_joule_integral > 0.0 &&
           summary.three_body_contact.conductor_joule_integral > 0.0 &&
           summary.conductor_joule_rel <= summary.spec.max_conductor_joule_rel &&
           summary.interface_spike.spike_hard_pass && summary.comparison_csv_written && summary.summary_csv_written;
}

inline void printTeam7SingleBodyVsThreeBodyContactComparisonSummary(
    const char *test_name, const AphiTeam7SingleBodyVsContactComparisonSummary &summary)
{
    std::cout << test_name << " single_body converged=" << (summary.single_body.converged ? 1 : 0)
              << " conductor_joule=" << summary.single_body.conductor_joule_integral
              << " three_body converged=" << (summary.three_body_contact.converged ? 1 : 0)
              << " conductor_joule=" << summary.three_body_contact.conductor_joule_integral
              << " conductor_joule_rel=" << summary.conductor_joule_rel
              << " b_probe_profile_rel=" << summary.b_probe_profile_rel
              << " b_probe_finite=" << (std::isfinite(summary.b_probe_profile_rel) ? 1 : 0)
              << " b_probe_informative_only=1"
              << " interface_spike_hard_pass=" << (summary.interface_spike.spike_hard_pass ? 1 : 0)
              << " comparison_csv_written=" << (summary.comparison_csv_written ? 1 : 0)
              << " summary_csv_written=" << (summary.summary_csv_written ? 1 : 0)
              << " report_dir=" << summary.spec.report_dir << " lambda_A=off" << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
