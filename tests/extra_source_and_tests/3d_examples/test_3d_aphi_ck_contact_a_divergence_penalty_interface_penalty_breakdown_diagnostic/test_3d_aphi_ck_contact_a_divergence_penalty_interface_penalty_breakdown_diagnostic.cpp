/**
 * Stage 10.12-B++: decompose Contact interface MMS defect into base operator vs A-penalty increment.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_interface_penalty_breakdown_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real x_interface = 0.5;
    const Real interface_band_half_width = 2.0 * dp_0;
    const Real max_interface_base_rel = 1.0e-3;
    const Real eta_values[] = {0.0, 0.1, 0.2};

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    AphiContactInterfacePenaltyBreakdownRow eta0_row;
    size_t eta_positive_reported = 0;

    for (const Real eta_a : eta_values)
    {
        const AphiContactInterfacePenaltyBreakdownRow row = runContactInterfacePenaltyBreakdownRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, scale_metrics, x_interface, interface_band_half_width);
        printContactInterfacePenaltyBreakdownRow(
            "test_3d_aphi_ck_contact_a_divergence_penalty_interface_penalty_breakdown_diagnostic", row);
        if (eta_a <= TinyReal)
        {
            eta0_row = row;
            continue;
        }
        if (std::isfinite(row.interface_penalty_only.relative_l2))
        {
            eta_positive_reported += 1;
        }
    }

    const bool eta0_ok = eta0_row.interface_base.relative_l2 <= max_interface_base_rel &&
                         eta0_row.interface_penalty_only.relative_l2 <= max_interface_base_rel;
    const bool passed = eta0_ok && eta_positive_reported == 2;

    std::cout << "test_3d_aphi_ck_contact_a_divergence_penalty_interface_penalty_breakdown_diagnostic eta0_ok="
              << (eta0_ok ? 1 : 0) << " eta_positive_reported=" << eta_positive_reported << "/2 passed="
              << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
