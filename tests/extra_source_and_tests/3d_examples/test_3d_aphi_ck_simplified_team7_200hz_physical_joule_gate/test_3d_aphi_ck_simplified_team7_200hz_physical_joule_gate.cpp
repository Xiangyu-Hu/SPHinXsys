/** Stage 10.16 pre-P1: 200 Hz physical Joule gate (Policy B, conductor-only physical heating). */
#include "electromagnetic_dynamics/diagnostics/aphi_simplified_team7_source_driven_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    static constexpr const char *k_test_name = "test_3d_aphi_ck_simplified_team7_200hz_physical_joule_gate";
    AphiSourceDrivenEmSolveSpec spec;
    spec.omega = 2.0 * Pi * 200.0;
    spec.write_vtp = false;
    spec.write_probe_csv = false;
    /**
     * Policy B: keep source-driven EM solve diagnostics unchanged and avoid
     * using air/conductor Joule ratio as a hard physical heating gate.
     */
    spec.max_air_to_conductor_joule_ratio = 1.0;

    const AphiSourceDrivenEmSolveMetrics metrics = runSourceDrivenEmSolve(ac, av, spec);
    const AphiPhysicalJoulePolicyBResult policy_b = evaluateTeam7PhysicalJoulePolicyB(metrics, spec);

    printSourceDrivenEmSolveMetrics(k_test_name, metrics, policy_b.physical_joule_gate_passed);
    std::cout << k_test_name
              << " policy=conductor_only_physical_joule"
              << " air_to_conductor_joule_ratio=" << policy_b.air_to_conductor_joule_ratio
              << " numerical_air_ratio_threshold=" << policy_b.numerical_air_ratio_threshold
              << " air_joule_numerical_artifact=" << (policy_b.air_joule_numerical_artifact ? 1 : 0)
              << " physical_heating_target=conductor_only"
              << " physical_joule_gate_passed=" << (policy_b.physical_joule_gate_passed ? 1 : 0) << std::endl;
    std::cout << k_test_name << " passed=" << (policy_b.physical_joule_gate_passed ? 1 : 0) << std::endl;
    return policy_b.physical_joule_gate_passed ? 0 : 1;
}
