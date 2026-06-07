/** Stage 10.16 P2: TEAM7 probe/reference pipeline alignment audit (diagnostic-only). */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_probe_reference_alignment_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_probe_reference_alignment_diagnostic";
    const AphiTeam7ProbeReferenceAlignmentSummary summary = runTeam7ProbeReferenceAlignmentAudit(ac, av);
    printTeam7ProbeReferenceAlignmentSummary(k_test_name, summary);
    const bool passed = team7ProbeReferenceAlignmentPassed(summary);
    std::cout << k_test_name << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
