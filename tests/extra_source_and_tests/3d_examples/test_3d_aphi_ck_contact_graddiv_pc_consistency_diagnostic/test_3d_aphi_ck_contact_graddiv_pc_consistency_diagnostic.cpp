/**
 * Stage 10.12-C2: split two-body Contact pairwise penalty apply FD vs JacobiGradDivABlock (Inner+Contact PC).
 */
#include "sphinxsys.h"
#include "electromagnetic_dynamics/diagnostics/aphi_contact_graddiv_pc_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real max_core_pc_diff = 5.0e-5;
    const AphiContactGradDivPcConsistencyMetrics metrics = runContactGradDivPcConsistencyMetrics(ac, av);
    printContactGradDivPcConsistencyMetrics("test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic", metrics);

    const bool passed = metrics.core_pc_probed > 0 && metrics.core_pc_max_abs_diff < max_core_pc_diff;

    std::cout << "test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
