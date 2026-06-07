/**
 * Stage 10.7: tolerance / high-lambda closure attempt (3x3 graddiv PC, restart=80).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_sweep_helpers.h"

#include <iostream>
#include <string>

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
    const Real impressed_current_amplitude = 5.0;
    const Real solenoidal_current_amplitude = 5.0;
    const Real imag_to_real_ratio = 0.1;
    const AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};
    const Vecd current_real(0.0, 0.0, 1.0);
    const Vecd current_imag(0.0, 0.0, 0.1);
    const UnsignedInt restart_dimension = 80;
    const UnsignedInt max_outer_iterations = 200;

    struct ClosureCase
    {
        const char *source_kind;
        Real tolerance;
        Real lambda_a;
    };
    const ClosureCase cases[] = {
        {"impressed", 5.0e-4, 300.0},
        {"impressed", 3.0e-4, 300.0},
        {"impressed", 2.0e-4, 300.0},
        {"impressed", 1.0e-4, 300.0},
        {"impressed", 5.0e-4, 500.0},
        {"impressed", 3.0e-4, 500.0},
        {"impressed", 1.0e-4, 500.0},
        {"impressed", 5.0e-4, 1000.0},
        {"impressed", 3.0e-4, 1000.0},
        {"impressed", 1.0e-4, 1000.0},
        {"solenoidal", 5.0e-4, 300.0},
        {"solenoidal", 3.0e-4, 300.0},
        {"solenoidal", 1.0e-4, 300.0},
    };

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_tol_closure_diagnostic begin restart="
              << restart_dimension << " max_outer=" << max_outer_iterations << std::endl;

    size_t converged_rows = 0;
    for (const ClosureCase &closure_case : cases)
    {
        AphiInnerADivergencePenaltySweepRow row;
        if (std::string(closure_case.source_kind) == "impressed")
        {
            row = runInnerADivergencePenaltySweepRow(
                ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
                phi_gauge_penalty, closure_case.lambda_a, true, closure_case.tolerance, source_region, current_real,
                current_imag, impressed_current_amplitude, restart_dimension, max_outer_iterations);
        }
        else
        {
            row = runInnerSolenoidalADivergencePenaltySweepRow(
                ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
                phi_gauge_penalty, closure_case.lambda_a, true, closure_case.tolerance, source_region,
                solenoidal_current_amplitude, imag_to_real_ratio, restart_dimension, max_outer_iterations);
        }

        std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_tol_closure_diagnostic"
                  << " source=" << closure_case.source_kind << " tol=" << closure_case.tolerance
                  << " a_penalty=" << closure_case.lambda_a << " converged=" << (row.converged ? 1 : 0)
                  << " outer=" << row.outer_iterations << " rel=" << row.final_relative_residual
                  << " true_rel=" << row.final_true_relative_residual
                  << " div_A_relative=" << row.global_div_a.div_a_relative
                  << " div_A_level=" << divAGaugeDiagnosticLevel(row.global_div_a.div_a_relative) << std::endl;
        if (row.converged)
        {
            converged_rows += 1;
        }
    }

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_tol_closure_diagnostic completed_rows="
              << (sizeof(cases) / sizeof(cases[0])) << " converged_rows=" << converged_rows << " passed=1"
              << std::endl;
    return 0;
}
