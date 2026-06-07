/**
 * Stage 10.7 inner: homogeneous impressed-current box, sweep lambda_A with post-solve divA.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_sweep_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real tolerance = 1.0e-4;
    const Real impressed_current_amplitude = 5.0;
    const AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};
    const Vecd current_real(0.0, 0.0, 1.0);
    const Vecd current_imag(0.0, 0.0, 0.1);

    struct SweepCase
    {
        bool use_a_divergence_penalty;
        Real a_divergence_penalty;
    };
    const SweepCase sweep_cases[] = {
        {false, 0.0},
        {true, 10.0},
        {true, 50.0},
        {true, 100.0},
        {true, 300.0},
        {true, 1000.0},
    };

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_sweep_diagnostic begin phi_gauge_penalty="
              << phi_gauge_penalty << std::endl;

    size_t converged_rows = 0;
    for (const SweepCase &sweep_case : sweep_cases)
    {
        const AphiInnerADivergencePenaltySweepRow row = runInnerADivergencePenaltySweepRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, sweep_case.a_divergence_penalty, sweep_case.use_a_divergence_penalty, tolerance,
            source_region, current_real, current_imag, impressed_current_amplitude);
        printInnerADivergencePenaltySweepRow("test_3d_aphi_ck_inner_a_divergence_penalty_sweep_diagnostic", row);
        if (row.converged)
        {
            converged_rows += 1;
        }
    }

    const bool passed = true;
    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_sweep_diagnostic completed_rows="
              << (sizeof(sweep_cases) / sizeof(sweep_cases[0])) << " converged_rows=" << converged_rows
              << " passed=" << (passed ? 1 : 0) << " (diagnostic: converged_rows not gating)" << std::endl;
    return 0;
}
