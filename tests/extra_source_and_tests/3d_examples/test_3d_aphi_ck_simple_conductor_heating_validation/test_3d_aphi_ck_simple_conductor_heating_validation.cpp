/**
 * Stage 10.14 P3: conductor heating energy check (uniform Joule, adiabatic lumped).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    AphiEmJouleThermalCouplingSpec spec;
    spec.thermal_mapping_only = true;
    spec.source_driven_em_thermal = false;
    const AphiEmJouleThermalCouplingMetrics metrics = runEmJouleThermalOneWayCoupling(ac, av, spec);

    const Real delta_T_avg = metrics.T_avg_end - metrics.T_avg_start;
    const Real body_volume = spec.body_length * spec.body_height * spec.body_width;
    const Real expected_delta_T_from_power = metrics.P_Joule * metrics.delta_t / (spec.rho_cp * body_volume + TinyReal);
    const Real delta_T_error = std::abs(delta_T_avg - expected_delta_T_from_power) /
                               (std::abs(expected_delta_T_from_power) + TinyReal);

    const bool passed = emJouleThermalCouplingPassed(metrics, spec) && delta_T_error < 0.2;

    std::cout << "test_3d_aphi_ck_simple_conductor_heating_validation"
              << " passed=" << (passed ? 1 : 0) << " Delta_T_avg=" << delta_T_avg
              << " expected_delta_T_avg=" << expected_delta_T_from_power << " delta_T_error=" << delta_T_error
              << " joule_source=controlled_uniform" << std::endl;
    return passed ? 0 : 1;
}
