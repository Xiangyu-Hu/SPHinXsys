/** Stage 10.18 entry: native TEAM7 eddy Joule post-process on SI reload (50 Hz, coil-path source). */
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_reload_geometry_helpers.h"

#include <filesystem>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_team7_native_joule_post_smoke";
    AphiTeam7NativeReloadSourceDrivenSmokeSpec spec;
    spec.omega = 2.0 * Pi * 50.0;
    spec.write_vtp = false;
    const AphiTeam7NativeReloadSourceDrivenSmokeMetrics metrics = runTeam7NativeReloadSourceDrivenSmoke(ac, av, spec);
    const bool passed = team7NativeReloadSourceDrivenSmokePassed(metrics, spec);

    std::error_code mkdir_error;
    std::filesystem::create_directories("team7_native_joule_post_report", mkdir_error);
    const bool csv_ok = writeTeam7NativeJoulePostSummaryCsv("team7_native_joule_post_report/summary.csv", metrics);

    printTeam7NativeReloadSourceDrivenSmoke(k_test_name, metrics, spec, passed);
    const Real air_to_conductor = metrics.air_joule_integral / (metrics.plate_joule_integral + TinyReal);
    std::cout << k_test_name << " joule_post_csv=" << (csv_ok ? 1 : 0)
              << " air_to_conductor_joule=" << air_to_conductor << " stage_10_18_entry=1" << std::endl;
    return (passed && csv_ok) ? 0 : 1;
}
