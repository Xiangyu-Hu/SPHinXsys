/**
 * Stage 10.13 P1/P2: cold-crucible-like A-phi electromagnetic demonstration case.
 *
 * Pipeline smoke test only — NOT quantitative validation.
 * MVP-A: analytic impressed vector potential A, then E/J/Joule/B/H post-processing + VTP.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_cold_crucible_demo_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    AphiColdCrucibleDemoSpec spec;
    const AphiColdCrucibleDemoMetrics metrics = runColdCrucibleImpressedDemo(ac, av, spec);
    const bool passed = coldCrucibleDemoMetricsPassed(metrics, spec);

    std::cout << "test_3d_aphi_ck_cold_crucible_demo"
              << " mode=impressed"
              << " dp=" << spec.dp << " particles=" << metrics.particles << " omega=" << spec.omega
              << " max_A=" << metrics.max_a << " max_B=" << metrics.max_b << " max_E=" << metrics.max_e
              << " max_J=" << metrics.max_j << " max_Joule=" << metrics.max_joule
              << " melt_Joule_integral=" << metrics.melt_joule_integral
              << " air_Joule_integral=" << metrics.air_joule_integral << " melt_J_max=" << metrics.melt_j_max
              << " air_J_max=" << metrics.air_j_max << " vtp_dir=" << metrics.vtp_output_dir
              << " demonstration_only=1 passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
