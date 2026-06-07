/** Stage 10.15: simplified TEAM7 source-driven dual-frequency (50 Hz / 200 Hz) diagnostic. */
#include "electromagnetic_dynamics/diagnostics/aphi_simplified_team7_source_driven_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_simplified_team7_dual_frequency_diagnostic";
    const AphiTeam7DualFrequencyComparisonSummary summary = runTeam7DualFrequencyComparison(ac, av);
    printTeam7DualFrequencyComparisonSummary(k_test_name, summary);
    const bool passed = team7DualFrequencyComparisonPassed(summary);
    std::cout << k_test_name << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
