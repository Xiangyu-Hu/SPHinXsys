/**
 * Stage 10.13-P4: Contact block-Jacobi graddiv PC with InnerOnly penalty stencil.
 * Verifies core FD vs PC (inner-only penalty pipeline) and that Contact graddiv is skipped on PC.
 */
#include "sphinxsys.h"
#include "electromagnetic_dynamics/diagnostics/aphi_contact_inner_only_graddiv_pc_consistency_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real max_core_pc_diff = 5.0e-5;
    const Real min_interface_pc_fd_gap = 1.0;
    const char *test_name = "test_3d_aphi_ck_contact_inner_only_graddiv_pc_consistency_diagnostic";

    const AphiContactInnerOnlyGradDivPcConsistencyResult result =
        runContactInnerOnlyGradDivPcConsistencyDiagnostic(ac, av, min_interface_pc_fd_gap);
    printContactInnerOnlyGradDivPcConsistencyResult(test_name, result);

    const bool passed = contactInnerOnlyGradDivPcConsistencyPassed(result, max_core_pc_diff);
    std::cout << test_name << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
