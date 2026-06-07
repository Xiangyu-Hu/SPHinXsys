/** Stage 10.16 P1: real annular prescribed-current source region (sigma=0, RHS-only). */
#include "electromagnetic_dynamics/diagnostics/aphi_real_annular_source_region_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_real_annular_source_region_diagnostic";
    const AphiRealAnnularSourceRegionSpec spec;
    const AphiSourceDrivenEmSolveMetrics metrics = runRealAnnularSourceRegionDiagnostic(ac, av, spec);
    const bool passed = realAnnularSourceRegionDiagnosticPassed(metrics, spec);
    printRealAnnularSourceRegionDiagnosticMetrics(k_test_name, metrics, spec, passed);
    std::cout << k_test_name << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
