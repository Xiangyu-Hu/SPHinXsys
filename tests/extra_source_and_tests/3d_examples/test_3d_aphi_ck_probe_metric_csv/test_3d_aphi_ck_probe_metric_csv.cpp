/** Stage 10.14 P7 (addendum) / P5 (main plan): TEAM7-oriented probe CSV infrastructure. */
#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"

using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    AphiSourceDrivenEmSolveSpec spec;
    spec.write_vtp = false;
    spec.write_probe_csv = true;
    const AphiSourceDrivenEmSolveMetrics metrics = runSourceDrivenEmSolve(ac, av, spec);

    const bool em_passed = sourceDrivenEmSolvePassed(metrics, spec);
    const bool probe_passed = metrics.probe_csv_requested && probeMetricCsvPassed(metrics.probe_csv);
    printSourceDrivenEmSolveMetrics("test_3d_aphi_ck_probe_metric_csv em", metrics, em_passed);
    printProbeCsvWriteResult("test_3d_aphi_ck_probe_metric_csv", metrics.probe_csv, probe_passed);

    const bool passed = em_passed && probe_passed;
    std::cout << "test_3d_aphi_ck_probe_metric_csv passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
