/**
 * Stage 10B: quantitative TEAM7-like reference comparison.
 * Reference = fine mesh (dp=0.075); candidate = engineering mesh (dp=0.1).
 * External FEM comparison deferred; uses self-reference + regression anchors.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

void printRunSummary(const char *label, const AphiTeam7CanonicalCaseRunResult &run)
{
    std::cout << " " << label << "_dp=" << run.dp << " particles=" << run.particles << " em_rel=" << run.metrics.em_rel
              << " outer=" << run.metrics.outer_iterations << " converged=" << (run.metrics.converged ? 1 : 0)
              << " conductor_joule=" << run.metrics.conductor_joule_power
              << " conductor_solution_block_max=" << run.metrics.conductor_solution_block_max
              << " coil_joule=" << run.metrics.coil_joule_power;
}

} // namespace

int main(int ac, char *av[])
{
    const AphiTeam7CanonicalCaseSpec spec;

    const AphiTeam7CanonicalCaseRunResult reference_run = runTeam7PhysicalCanonicalCase(ac, av, spec.reference_dp);
    const AphiTeam7CanonicalCaseRunResult candidate_run = runTeam7PhysicalCanonicalCase(ac, av, spec.candidate_dp);

    const Real conductor_joule_rel =
        relativeMetricChange(reference_run.metrics.conductor_joule_power, candidate_run.metrics.conductor_joule_power);
    const Real conductor_solution_block_rel =
        relativeMetricChange(reference_run.metrics.conductor_solution_block_max,
                             candidate_run.metrics.conductor_solution_block_max);
    const Real centerline_profile_rel =
        hostCenterlineProfileL2RelativeDifference(reference_run.centerline, candidate_run.centerline);

    const Real reference_joule_recorded_rel =
        relativeMetricChange(spec.recorded_conductor_joule_power, reference_run.metrics.conductor_joule_power);
    const Real reference_solution_block_recorded_rel =
        relativeMetricChange(spec.recorded_conductor_solution_block_max,
                             reference_run.metrics.conductor_solution_block_max);

    const bool reference_converged = reference_run.metrics.converged;
    const bool candidate_converged = candidate_run.metrics.converged;
    const bool candidate_vs_reference_passed =
        conductor_joule_rel <= spec.max_candidate_vs_reference_rel &&
        conductor_solution_block_rel <= spec.max_candidate_vs_reference_rel;
    const bool centerline_informative_ok = centerline_profile_rel <= spec.max_centerline_profile_rel;
    const bool reference_regression_passed =
        reference_joule_recorded_rel <= spec.max_reference_vs_recorded_rel &&
        reference_solution_block_recorded_rel <= spec.max_reference_vs_recorded_rel;

    const bool passed = reference_converged && candidate_converged && candidate_vs_reference_passed &&
                        reference_regression_passed;

    std::cout << "test_3d_aphi_ck_gmres_team7_quantitative_reference_comparison";
    printRunSummary("reference", reference_run);
    printRunSummary("candidate", candidate_run);
    std::cout << " conductor_joule_rel=" << conductor_joule_rel
              << " conductor_solution_block_rel=" << conductor_solution_block_rel
              << " centerline_profile_rel=" << centerline_profile_rel
              << " reference_joule_recorded_rel=" << reference_joule_recorded_rel
              << " reference_solution_block_recorded_rel=" << reference_solution_block_recorded_rel
              << " candidate_vs_reference_passed=" << (candidate_vs_reference_passed ? 1 : 0)
              << " centerline_informative_ok=" << (centerline_informative_ok ? 1 : 0)
              << " reference_regression_passed=" << (reference_regression_passed ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
