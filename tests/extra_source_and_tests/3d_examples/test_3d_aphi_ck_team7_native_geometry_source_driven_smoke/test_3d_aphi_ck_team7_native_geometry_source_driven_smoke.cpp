/** Stage 10.17: source-driven A-phi smoke on native TEAM7 reload particles (SYCL CK). */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_reload_geometry_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_native_geometry_source_driven_smoke";
    AphiTeam7NativeReloadSourceDrivenSmokeSpec spec;
    spec.omega = 2.0 * Pi * 50.0;
    const AphiTeam7NativeReloadSourceDrivenSmokeMetrics metrics = runTeam7NativeReloadSourceDrivenSmoke(ac, av, spec);
    const bool passed = team7NativeReloadSourceDrivenSmokePassed(metrics, spec);
    printTeam7NativeReloadSourceDrivenSmoke(k_test_name, metrics, spec, passed);
    std::cout << k_test_name << " native_reload_geometry=1 standard_team7_validation=0" << std::endl;
    return passed ? 0 : 1;
}
