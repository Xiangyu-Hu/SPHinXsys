/**
 * Stage 10.12-B+: spatial breakdown of Contact two-body MMS consistency defect ||Ku-b||/||b||.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_interface_mms_consistency_helpers.h"

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
    const Real max_global_rel_l2_eta0 = 1.0e-4;
    const Real max_core_safe_rel_l2_eta0 = 1.0e-4;
    const Real max_interface_rel_l2_eta_positive = 0.15;
    const Real eta_values[] = {0.0, 0.1, 0.2};

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    AphiContactInterfaceMmsConsistencyRow eta0_row;
    size_t eta_positive_reported = 0;

    for (const Real eta_a : eta_values)
    {
        const AphiContactInterfaceMmsConsistencyRow row = runContactInterfaceMmsConsistencyRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, eta_a, scale_metrics, x_interface, interface_band_half_width,
            AphiDivFreeValidationFieldKind::Az2D);
        printContactInterfaceMmsConsistencyRow(
            "test_3d_aphi_ck_contact_a_divergence_penalty_interface_mms_consistency_diagnostic", row);
        if (eta_a <= TinyReal)
        {
            eta0_row = row;
            continue;
        }
        if (std::isfinite(row.global_relative_l2) && row.interface_band.particle_count > 0)
        {
            eta_positive_reported += 1;
        }
    }

    const bool eta0_ok = eta0_row.global_relative_l2 <= max_global_rel_l2_eta0 &&
                         eta0_row.stencil_safe_core.relative_l2 <= max_core_safe_rel_l2_eta0;
    const bool passed = eta0_ok && eta_positive_reported == 2;

    std::cout << "test_3d_aphi_ck_contact_a_divergence_penalty_interface_mms_consistency_diagnostic eta0_ok="
              << (eta0_ok ? 1 : 0) << " eta_positive_reported=" << eta_positive_reported << "/2 passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
