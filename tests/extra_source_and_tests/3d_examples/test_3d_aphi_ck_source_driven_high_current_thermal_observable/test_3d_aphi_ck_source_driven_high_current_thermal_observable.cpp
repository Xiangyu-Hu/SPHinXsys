/** Stage 10.16 P5: source-driven EM Joule -> thermal with elevated impressed current (no uniform Joule). */
#include "electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_source_driven_high_current_thermal_observable";

    AphiEmJouleThermalCouplingSpec spec;
    spec.thermal_use_uniform_joule = false;
    spec.source_driven_em_thermal = true;
    spec.impressed_current_amplitude = 40.0 * benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    spec.source_driven_delta_t = 0.05;
    spec.source_driven_thermal_steps = 80;
    spec.max_source_driven_energy_relative_error = 0.08;
    spec.write_vtp = false;

    const AphiEmJouleThermalCouplingMetrics metrics = runEmJouleThermalOneWayCoupling(ac, av, spec);
    const bool passed = sourceDrivenHighCurrentThermalObservablePassed(metrics, spec);

    printEmJouleThermalCouplingMetrics(k_test_name, metrics, passed);
    std::cout << k_test_name << " impressed_current_amplitude=" << spec.impressed_current_amplitude
              << " thermal_use_uniform_joule=0 source_driven_em_thermal=1"
              << " high_current_thermal_observable_passed=" << (passed ? 1 : 0) << std::endl;
    std::cout << k_test_name << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
