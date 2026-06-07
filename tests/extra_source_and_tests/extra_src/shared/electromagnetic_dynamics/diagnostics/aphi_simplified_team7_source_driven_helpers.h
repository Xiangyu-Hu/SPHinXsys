#ifndef APHI_SIMPLIFIED_TEAM7_SOURCE_DRIVEN_HELPERS_H
#define APHI_SIMPLIFIED_TEAM7_SOURCE_DRIVEN_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiTeam7DualFrequencyRunResult
{
    Real frequency_hz = 0.0;
    Real omega = 0.0;
    AphiSourceDrivenEmSolveSpec spec{};
    AphiSourceDrivenEmSolveMetrics metrics{};
    bool em_passed = false;
};

struct AphiTeam7DualFrequencyComparisonSummary
{
    std::string report_dir = "team7_dual_frequency_report";
    AphiTeam7DualFrequencyRunResult run_50hz{};
    AphiTeam7DualFrequencyRunResult run_200hz{};
    AphiProbeCsvWriteResult probe_csv_50hz{};
    AphiProbeCsvWriteResult probe_csv_200hz{};
    Real conductor_joule_ratio_200_to_50 = 0.0;
    Real max_abs_E_ratio_200_to_50 = 0.0;
    Real max_abs_J_conductor_ratio_200_to_50 = 0.0;
    bool frequency_response_informative = false;
    bool comparison_csv_written = false;
    bool summary_csv_written = false;
    bool probe_csv_written = false;
};

struct AphiPhysicalJoulePolicyBResult
{
    Real air_to_conductor_joule_ratio = 0.0;
    Real numerical_air_ratio_threshold = 0.05;
    bool air_joule_numerical_artifact = false;
    bool physical_heating_target_is_conductor_only = true;
    bool physical_joule_gate_passed = false;
};

inline AphiSourceDrivenEmSolveMetrics runSimplifiedTeam7SourceDriven(int ac, char *av[])
{
    AphiSourceDrivenEmSolveSpec spec;
    spec.omega = 2.0 * Pi * 50.0;
    spec.write_vtp = true;
    return runSourceDrivenEmSolve(ac, av, spec);
}

inline bool simplifiedTeam7SourceDrivenPassed(const AphiSourceDrivenEmSolveMetrics &metrics,
                                              const AphiSourceDrivenEmSolveSpec &spec)
{
    return sourceDrivenEmSolvePassed(metrics, spec);
}

inline Real airToConductorJouleRatio(const AphiSourceDrivenEmSolveMetrics &metrics)
{
    return metrics.air_Joule_integral / (metrics.conductor_Joule_integral + TinyReal);
}

inline bool writeTeam7DualFrequencyComparisonCsv(const std::string &path,
                                                 const AphiTeam7DualFrequencyRunResult &run_50hz,
                                                 const AphiTeam7DualFrequencyRunResult &run_200hz)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "frequency_hz,omega,em_passed,converged,num_iterations,final_residual,max_abs_A,max_abs_E,"
              "max_abs_J_conductor,max_Joule_conductor,conductor_Joule,air_Joule,source_Joule,"
              "air_to_conductor_joule,source_rhs_l2_source_region,probe_output_dir\n";

    const auto write_row = [&](const AphiTeam7DualFrequencyRunResult &run) {
        const AphiSourceDrivenEmSolveMetrics &m = run.metrics;
        output << run.frequency_hz << "," << run.omega << "," << (run.em_passed ? 1 : 0) << ","
               << (m.converged ? 1 : 0) << "," << m.num_iterations << "," << m.final_residual << "," << m.max_abs_A
               << "," << m.max_abs_E << "," << m.max_abs_J_conductor << "," << m.max_Joule_conductor << ","
               << m.conductor_Joule_integral << "," << m.air_Joule_integral << "," << m.source_Joule_integral << ","
               << airToConductorJouleRatio(m) << "," << m.source_rhs_l2_source_region << ","
               << run.spec.probe_output_dir << "\n";
    };
    write_row(run_50hz);
    write_row(run_200hz);
    return true;
}

inline bool writeTeam7DualFrequencySummaryCsv(const std::string &path,
                                              const AphiTeam7DualFrequencyComparisonSummary &summary)
{
    std::ofstream output(path);
    if (!output)
    {
        return false;
    }
    output << std::setprecision(10);
    output << "metric,value\n";
    output << "conductor_joule_ratio_200_to_50," << summary.conductor_joule_ratio_200_to_50 << "\n";
    output << "max_abs_E_ratio_200_to_50," << summary.max_abs_E_ratio_200_to_50 << "\n";
    output << "max_abs_J_conductor_ratio_200_to_50," << summary.max_abs_J_conductor_ratio_200_to_50 << "\n";
    output << "frequency_response_informative," << (summary.frequency_response_informative ? 1 : 0) << "\n";
    output << "50hz_em_passed," << (summary.run_50hz.em_passed ? 1 : 0) << "\n";
    output << "200hz_em_passed," << (summary.run_200hz.em_passed ? 1 : 0) << "\n";
    return true;
}

inline AphiTeam7DualFrequencyRunResult runTeam7DualFrequencySourceDrivenCase(int ac, char *av[], Real frequency_hz,
                                                                             Real max_air_to_conductor_joule_ratio,
                                                                             const std::string &probe_output_dir,
                                                                             bool write_probe_csv)
{
    AphiTeam7DualFrequencyRunResult result;
    result.frequency_hz = frequency_hz;
    result.omega = 2.0 * Pi * frequency_hz;
    result.spec.omega = result.omega;
    result.spec.write_vtp = false;
    result.spec.write_probe_csv = write_probe_csv;
    result.spec.probe_output_dir = probe_output_dir;
    result.spec.max_air_to_conductor_joule_ratio = max_air_to_conductor_joule_ratio;
    if (write_probe_csv)
    {
        std::error_code mkdir_error;
        std::filesystem::create_directories(probe_output_dir, mkdir_error);
    }
    result.metrics = runSourceDrivenEmSolve(ac, av, result.spec);
    result.em_passed = sourceDrivenEmSolvePassed(result.metrics, result.spec);
    return result;
}

inline AphiTeam7DualFrequencyComparisonSummary runTeam7DualFrequencyComparison(
    int ac, char *av[], const std::string &report_dir = "team7_dual_frequency_report", bool write_probe_csv = true)
{
    std::error_code mkdir_error;
    std::filesystem::create_directories(report_dir, mkdir_error);

    AphiTeam7DualFrequencyComparisonSummary summary;
    summary.report_dir = report_dir;
    // Keep probe/report artifacts outside SPH's mutable `output/` folder (relocated to output_backup each run).
    summary.run_50hz = runTeam7DualFrequencySourceDrivenCase(ac, av, 50.0, 0.05, report_dir + "/probe_50hz",
                                                             write_probe_csv);
    summary.run_200hz = runTeam7DualFrequencySourceDrivenCase(ac, av, 200.0, 1.0, report_dir + "/probe_200hz",
                                                                write_probe_csv);
    summary.conductor_joule_ratio_200_to_50 =
        summary.run_200hz.metrics.conductor_Joule_integral /
        (summary.run_50hz.metrics.conductor_Joule_integral + TinyReal);
    summary.max_abs_E_ratio_200_to_50 = summary.run_200hz.metrics.max_abs_E /
                                        (summary.run_50hz.metrics.max_abs_E + TinyReal);
    summary.max_abs_J_conductor_ratio_200_to_50 = summary.run_200hz.metrics.max_abs_J_conductor /
                                                  (summary.run_50hz.metrics.max_abs_J_conductor + TinyReal);
    summary.frequency_response_informative =
        std::isfinite(summary.conductor_joule_ratio_200_to_50) && summary.conductor_joule_ratio_200_to_50 > 0.0;
    summary.probe_csv_50hz = summary.run_50hz.metrics.probe_csv;
    summary.probe_csv_200hz = summary.run_200hz.metrics.probe_csv;
    summary.comparison_csv_written = writeTeam7DualFrequencyComparisonCsv(
        report_dir + "/aphi_team7_dual_frequency_comparison.csv", summary.run_50hz, summary.run_200hz);
    summary.summary_csv_written =
        writeTeam7DualFrequencySummaryCsv(report_dir + "/aphi_team7_dual_frequency_summary.csv", summary);
    summary.probe_csv_written =
        !write_probe_csv ||
        (probeMetricCsvPassed(summary.probe_csv_50hz) && probeMetricCsvPassed(summary.probe_csv_200hz));
    return summary;
}

inline bool team7DualFrequencyComparisonPassed(const AphiTeam7DualFrequencyComparisonSummary &summary)
{
    return summary.run_50hz.em_passed && summary.run_200hz.em_passed && summary.frequency_response_informative &&
           summary.comparison_csv_written && summary.summary_csv_written && summary.probe_csv_written;
}

inline AphiPhysicalJoulePolicyBResult evaluateTeam7PhysicalJoulePolicyB(const AphiSourceDrivenEmSolveMetrics &metrics,
                                                                         const AphiSourceDrivenEmSolveSpec &spec)
{
    AphiPhysicalJoulePolicyBResult result;
    result.air_to_conductor_joule_ratio = airToConductorJouleRatio(metrics);
    result.air_joule_numerical_artifact = result.air_to_conductor_joule_ratio > result.numerical_air_ratio_threshold;
    result.physical_heating_target_is_conductor_only = true;
    result.physical_joule_gate_passed = metrics.converged && metrics.finite_field_check &&
                                        metrics.max_abs_A > spec.min_solution_norm &&
                                        metrics.max_abs_E > spec.min_conductor_E &&
                                        metrics.max_abs_J_conductor > spec.min_conductor_J &&
                                        metrics.max_Joule_conductor > 0.0 &&
                                        metrics.conductor_Joule_integral > spec.min_conductor_joule_integral &&
                                        metrics.source_rhs_l2_source_region > 0.0 &&
                                        metrics.particle_count_source > 0 && metrics.particle_count_conductor > 0 &&
                                        metrics.particle_count_air > 0;
    return result;
}

inline void printTeam7DualFrequencyComparisonSummary(const char *test_name,
                                                   const AphiTeam7DualFrequencyComparisonSummary &summary)
{
    const std::string label_50hz = std::string(test_name) + " 50hz";
    const std::string label_200hz = std::string(test_name) + " 200hz";
    printSourceDrivenEmSolveMetrics(label_50hz.c_str(), summary.run_50hz.metrics, summary.run_50hz.em_passed);
    printSourceDrivenEmSolveMetrics(label_200hz.c_str(), summary.run_200hz.metrics, summary.run_200hz.em_passed);
    std::cout << test_name << " conductor_joule_ratio_200_to_50=" << summary.conductor_joule_ratio_200_to_50
              << " max_abs_E_ratio_200_to_50=" << summary.max_abs_E_ratio_200_to_50
              << " max_abs_J_conductor_ratio_200_to_50=" << summary.max_abs_J_conductor_ratio_200_to_50
              << " frequency_response_informative=" << (summary.frequency_response_informative ? 1 : 0)
              << " comparison_csv_written=" << (summary.comparison_csv_written ? 1 : 0)
              << " summary_csv_written=" << (summary.summary_csv_written ? 1 : 0)
              << " probe_csv_written=" << (summary.probe_csv_written ? 1 : 0)
              << " diagnostic_only=1"
              << " report_dir=" << summary.report_dir << " probe_csv_50hz_dir=" << summary.run_50hz.spec.probe_output_dir
              << " probe_csv_200hz_dir=" << summary.run_200hz.spec.probe_output_dir
              << " simplified_team7=1 standard_team7_validation=0" << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
