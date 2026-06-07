/**
 * Stage 10.14 P1: source-driven A-phi GMRES (impressed coil RHS), E/J/Joule postprocess + VTP.
 * Not impressed-A assignment; lambda_A off.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const AphiSourceDrivenEmSolveSpec spec;
    const AphiSourceDrivenEmSolveMetrics metrics = runSourceDrivenEmSolve(ac, av, spec);
    const bool passed = sourceDrivenEmSolvePassed(metrics, spec);
    printSourceDrivenEmSolveMetrics("test_3d_aphi_ck_source_driven_em_solve", metrics, passed);
    return passed ? 0 : 1;
}
