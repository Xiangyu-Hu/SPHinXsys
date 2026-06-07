/**
 * Stage 10.14 P2: one-way Joule -> adiabatic lumped thermal (mapping-only + source-driven EM).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    AphiEmJouleThermalCouplingSpec mapping_spec;
    mapping_spec.thermal_mapping_only = true;
    mapping_spec.source_driven_em_thermal = false;
    const AphiEmJouleThermalCouplingMetrics mapping_metrics = runEmJouleThermalOneWayCoupling(ac, av, mapping_spec);
    const bool mapping_passed = emJouleThermalCouplingPassed(mapping_metrics, mapping_spec);
    printEmJouleThermalCouplingMetrics("test_3d_aphi_ck_em_joule_thermal_one_way mapping", mapping_metrics,
                                       mapping_passed);

    AphiEmJouleThermalCouplingSpec low_power_spec;
    low_power_spec.thermal_mapping_only = false;
    low_power_spec.source_driven_em_thermal = true;
    const AphiEmJouleThermalCouplingMetrics low_power_metrics = runEmJouleThermalOneWayCoupling(ac, av, low_power_spec);
    const bool low_power_gate_passed = emJouleThermalCouplingPassed(low_power_metrics, low_power_spec);
    const bool low_power_passed = low_power_gate_passed && low_power_metrics.resolution_floor_triggered;
    printEmJouleThermalCouplingMetrics("test_3d_aphi_ck_em_joule_thermal_one_way source_driven_low_power",
                                       low_power_metrics, low_power_passed);

    AphiEmJouleThermalCouplingSpec high_power_spec;
    high_power_spec.thermal_use_uniform_joule = true;
    high_power_spec.source_driven_em_thermal = false;
    high_power_spec.uniform_joule_source = 1.0e5;
    high_power_spec.delta_t = 0.01;
    high_power_spec.thermal_steps = 10;
    const AphiEmJouleThermalCouplingMetrics high_power_metrics =
        runEmJouleThermalOneWayCoupling(ac, av, high_power_spec);
    const bool high_power_gate_passed = emJouleThermalCouplingPassed(high_power_metrics, high_power_spec);
    const bool high_power_observability_passed =
        !high_power_metrics.resolution_floor_triggered &&
        high_power_metrics.observed_temperature_delta_max > high_power_metrics.temperature_resolution_floor;
    const bool high_power_passed = high_power_gate_passed && high_power_observability_passed;
    printEmJouleThermalCouplingMetrics("test_3d_aphi_ck_em_joule_thermal_one_way source_driven_high_power",
                                       high_power_metrics, high_power_passed);
    std::cout << "test_3d_aphi_ck_em_joule_thermal_one_way source_driven_high_power_observability"
              << " expected_no_floor=1"
              << " passed=" << (high_power_observability_passed ? 1 : 0) << std::endl;

    const bool passed = mapping_passed && low_power_passed && high_power_passed;
    std::cout << "test_3d_aphi_ck_em_joule_thermal_one_way passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
