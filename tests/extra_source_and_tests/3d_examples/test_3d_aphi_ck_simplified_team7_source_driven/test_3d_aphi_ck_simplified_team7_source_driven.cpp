/** Stage 10.14 P6: simplified TEAM7 source-driven (50 Hz), not standard TEAM7 validation. */
#include "electromagnetic_dynamics/diagnostics/aphi_simplified_team7_source_driven_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    AphiSourceDrivenEmSolveSpec spec;
    spec.omega = 2.0 * Pi * 50.0;
    const AphiSourceDrivenEmSolveMetrics metrics = runSourceDrivenEmSolve(ac, av, spec);
    const bool passed = sourceDrivenEmSolvePassed(metrics, spec);
    printSourceDrivenEmSolveMetrics("test_3d_aphi_ck_simplified_team7_source_driven", metrics, passed);
    std::cout << "test_3d_aphi_ck_simplified_team7_source_driven simplified_team7=1 standard_team7_validation=0"
              << std::endl;
    return passed ? 0 : 1;
}
