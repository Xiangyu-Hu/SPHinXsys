/** Stage 10.17: native TEAM7 reload geometry audit (no A-phi solve). */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_reload_geometry_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_native_geometry_audit";
    const AphiTeam7NativeGeometryAuditSummary summary = runTeam7NativeGeometryAudit(ac, av);
    printTeam7NativeGeometryAudit(k_test_name, summary);
    const bool passed = team7NativeGeometryAuditPassed(summary);
    std::cout << k_test_name << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
