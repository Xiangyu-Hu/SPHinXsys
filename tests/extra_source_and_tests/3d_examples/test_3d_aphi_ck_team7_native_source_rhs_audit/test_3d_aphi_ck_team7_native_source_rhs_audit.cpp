/** Stage 10.17 P1: TEAM7 native coil path RHS definition audit (SI reload, NI=2742 A). */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_reload_geometry_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_native_source_rhs_audit";
    Team7CoilPathSourceSpec coil_spec;
    const AphiTeam7NativeSourceRhsAuditSummary summary = runTeam7NativeSourceRhsAudit(ac, av, coil_spec);
    const bool passed = team7NativeSourceRhsAuditPassed(summary);
    printTeam7NativeSourceRhsAudit(k_test_name, summary, passed);
    return passed ? 0 : 1;
}
