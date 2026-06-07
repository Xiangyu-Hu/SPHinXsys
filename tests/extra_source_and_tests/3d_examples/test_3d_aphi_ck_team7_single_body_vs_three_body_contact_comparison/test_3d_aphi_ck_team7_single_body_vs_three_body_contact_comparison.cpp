/** Stage 10.16 P4: single-body region-tagged TEAM7 vs three-body Contact comparison. */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_single_body_vs_contact_comparison_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_single_body_vs_three_body_contact_comparison";
    const AphiTeam7SingleBodyVsContactComparisonSummary summary = runTeam7SingleBodyVsThreeBodyContactComparison(ac, av);
    printTeam7SingleBodyVsThreeBodyContactComparisonSummary(k_test_name, summary);
    const bool passed = team7SingleBodyVsThreeBodyContactComparisonPassed(summary);
    std::cout << k_test_name << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
