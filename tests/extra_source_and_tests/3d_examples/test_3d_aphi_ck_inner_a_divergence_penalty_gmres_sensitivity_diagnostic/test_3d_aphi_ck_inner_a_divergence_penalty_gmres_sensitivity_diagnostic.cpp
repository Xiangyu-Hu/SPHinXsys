/**
 * Stage 10.7: GMRES outer/restart sensitivity for inner A-divergence penalty (3x3 graddiv PC).
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
    const Real solenoidal_current_amplitude = 5.0;
    const Real imag_to_real_ratio = 0.1;
    const AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};
    const Vecd current_real(0.0, 0.0, 1.0);
    const Vecd current_imag(0.0, 0.0, 0.1);

    struct SolverSpec
    {
        UnsignedInt restart_dimension;
        UnsignedInt max_outer_iterations;
    };
    const SolverSpec solver_specs[] = {
        {50, 100},
        {50, 200},
        {50, 300},
        {80, 200},
    };
    const Real lambda_values[] = {50.0, 100.0, 300.0};

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic begin tolerance="
              << tolerance << " phi_gauge_penalty=" << phi_gauge_penalty << std::endl;

    size_t converged_rows = 0;
    size_t total_rows = 0;

    for (const Real lambda_a : lambda_values)
    {
        for (const SolverSpec &solver_spec : solver_specs)
        {
            const AphiInnerADivergencePenaltySweepRow impressed_row = runInnerADivergencePenaltySweepRow(
                ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
                phi_gauge_penalty, lambda_a, true, tolerance, source_region, current_real, current_imag,
                impressed_current_amplitude, solver_spec.restart_dimension, solver_spec.max_outer_iterations);
            printInnerADivergencePenaltySweepRow(
                "test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic impressed", impressed_row);
            total_rows += 1;
            if (impressed_row.converged)
            {
                converged_rows += 1;
            }
        }
    }

    for (const Real lambda_a : lambda_values)
    {
        for (const SolverSpec &solver_spec : solver_specs)
        {
            const AphiInnerADivergencePenaltySweepRow solenoidal_row = runInnerSolenoidalADivergencePenaltySweepRow(
                ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
                phi_gauge_penalty, lambda_a, true, tolerance, source_region, solenoidal_current_amplitude,
                imag_to_real_ratio, solver_spec.restart_dimension, solver_spec.max_outer_iterations);
            printInnerADivergencePenaltySweepRow(
                "test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic solenoidal", solenoidal_row);
            total_rows += 1;
            if (solenoidal_row.converged)
            {
                converged_rows += 1;
            }
        }
    }

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_gmres_sensitivity_diagnostic completed_rows=" << total_rows
              << " converged_rows=" << converged_rows << " passed=1 (diagnostic: converged_rows not gating)"
              << std::endl;
    return 0;
}
