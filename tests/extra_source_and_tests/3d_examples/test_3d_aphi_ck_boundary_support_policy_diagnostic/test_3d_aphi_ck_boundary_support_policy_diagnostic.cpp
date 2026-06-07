/** Stage 10.14 P5: outer padding sensitivity (fixed TEAM7 box, enlarged boundary_width). */
#include "electromagnetic_dynamics/diagnostics/aphi_boundary_support_policy_diagnostic_helpers.h"

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    StdVec<AphiBoundaryPolicyRow> rows;
    rows.push_back(runBoundaryPolicyRow(ac, av, AphiBoundarySupportPolicy::Baseline, 1.0));
    rows.push_back(runBoundaryPolicyRow(ac, av, AphiBoundarySupportPolicy::EnlargedAir, 2.0));
    rows.push_back(runBoundaryPolicyRow(ac, av, AphiBoundarySupportPolicy::EnlargedAir, 3.0));
    rows.push_back(runBoundaryPolicyRow(ac, av, AphiBoundarySupportPolicy::PassiveAirShell, 1.0, 2.0));
    const bool passed = boundarySupportPolicyDiagnosticPassed(rows);
    printBoundaryPolicyRows("test_3d_aphi_ck_boundary_support_policy_diagnostic", rows, passed);
    return passed ? 0 : 1;
}
