/** Stage 10.14-B P6: annular/racetrack source region smoke (RHS-only coil sigma=0). */
#include "electromagnetic_dynamics/diagnostics/aphi_annular_source_region_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const AphiAnnularSourceRegionSpec spec;
    const AphiSourceDrivenEmSolveMetrics metrics = runAnnularSourceRegionDiagnostic(ac, av, spec);
    const bool passed = annularSourceRegionDiagnosticPassed(metrics);
    printAnnularSourceRegionDiagnosticMetrics("test_3d_aphi_ck_annular_source_region_diagnostic", metrics, passed);
    return passed ? 0 : 1;
}
