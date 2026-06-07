/**
 * Stage 10C: cold-crucible geometry scaffold — EM-only solve + Joule source (no thermal feedback).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const AphiColdCrucibleCaseSpec spec;
    const AphiColdCrucibleScaffoldRunResult run = runColdCrucibleScaffoldCase(ac, av, spec.dp);

    const bool region_tagging_ok = run.region_counts.melt > 0 && run.region_counts.crucible_wall > 0 &&
                                   run.region_counts.coil > 0 && run.region_counts.air > 0;
    const bool passed = run.metrics.converged && region_tagging_ok && run.metrics.melt_joule_power > 0.0 &&
                        run.metrics.melt_solution_block_max >= spec.min_melt_solution_block_max &&
                        run.metrics.crucible_joule_power >= 0.0 && run.metrics.min_joule_source >= -1.0e-12;

    std::cout << "test_3d_aphi_ck_gmres_cold_crucible_scaffold"
              << " dp=" << run.dp << " particles=" << run.particles << " em_rel=" << run.metrics.em_rel
              << " outer=" << run.metrics.outer_iterations << " converged=" << (run.metrics.converged ? 1 : 0)
              << " region_melt=" << run.region_counts.melt << " region_crucible=" << run.region_counts.crucible_wall
              << " region_coil=" << run.region_counts.coil << " region_air=" << run.region_counts.air
              << " melt_joule=" << run.metrics.melt_joule_power
              << " crucible_joule=" << run.metrics.crucible_joule_power << " coil_joule=" << run.metrics.coil_joule_power
              << " total_joule=" << run.metrics.total_joule_power
              << " melt_solution_block_max=" << run.metrics.melt_solution_block_max
              << " crucible_solution_block_max=" << run.metrics.crucible_solution_block_max
              << " coil_solution_block_max=" << run.metrics.coil_solution_block_max
              << " min_joule=" << run.metrics.min_joule_source << " region_tagging_ok=" << (region_tagging_ok ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
