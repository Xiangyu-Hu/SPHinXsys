/**
 * Stage 10.11-A: boundary-layer divA diagnostic with core/boundary gradDen split + boundary_width sweep.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_boundary_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real dp_0 = 0.1;
    const Real core_shell = 2.5 * dp_0;
    const Real max_core_pairwise_grad_den = 1.0e-5;
    const Real max_core_b_grad_den = 0.05;
    const Real boundary_width_scales[] = {2.0, 3.0, 4.0};
    const AphiDivFreeValidationFieldKind field_kinds[] = {
        AphiDivFreeValidationFieldKind::Az2D,
        AphiDivFreeValidationFieldKind::CrossSine3D,
        AphiDivFreeValidationFieldKind::Sinusoidal3DCurlPsi,
    };

    size_t core_pass = 0;
    size_t core_expected = 0;
    AphiPairwiseDivABoundaryMetrics az2d_ref{};

    for (const Real width_scale : boundary_width_scales)
    {
        const Real boundary_width = width_scale * dp_0;
        for (const AphiDivFreeValidationFieldKind field_kind : field_kinds)
        {
            const AphiPairwiseDivABoundaryMetrics metrics = runPairwiseDivABoundaryMetrics(
                ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, field_kind);
            printPairwiseDivABoundaryMetrics("test_3d_aphi_ck_pairwise_diva_boundary_diagnostic", metrics);

            if (std::abs(width_scale - 3.0) > 1.0e-6)
            {
                continue;
            }
            if (field_kind == AphiDivFreeValidationFieldKind::Az2D)
            {
                az2d_ref = metrics;
            }
            if (field_kind == AphiDivFreeValidationFieldKind::Az2D ||
                field_kind == AphiDivFreeValidationFieldKind::CrossSine3D)
            {
                core_expected += 1;
                if (metrics.core_pairwise_grad_den.div_a_relative <= max_core_pairwise_grad_den)
                {
                    core_pass += 1;
                }
            }
        }
    }

    const bool az2d_boundary_dominates = az2d_ref.boundary_div_energy_fraction > 0.9;
    const bool az2d_core_b_ok = az2d_ref.core_b_corrected_grad_den.div_a_relative <= max_core_b_grad_den;
    const bool passed = core_pass == core_expected && az2d_boundary_dominates && az2d_core_b_ok;

    std::cout << "test_3d_aphi_ck_pairwise_diva_boundary_diagnostic core_pass=" << core_pass << "/"
              << core_expected << " az2d_boundary_energy_frac=" << az2d_ref.boundary_div_energy_fraction
              << " az2d_core_B_gradDen_rel=" << az2d_ref.core_b_corrected_grad_den.div_a_relative
              << " az2d_boundary_B_gradDen_rel=" << az2d_ref.boundary_b_corrected_grad_den.div_a_relative
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
