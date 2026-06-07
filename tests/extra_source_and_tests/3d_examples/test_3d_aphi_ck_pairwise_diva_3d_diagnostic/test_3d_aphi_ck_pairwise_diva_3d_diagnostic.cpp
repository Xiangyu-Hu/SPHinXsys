/**
 * Stage 10.10-B: pairwise vs B-corrected divA on four divergence-free validation fields + dp sweep.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h"

#include <cmath>
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
    const Real core_shell_scale = 2.5;
    const Real max_pairwise_rel_linear2d = 0.05;
    const Real max_pairwise_rel_az2d = 0.15;
    const Real max_pairwise_rel_crosssine3d = 0.15;
    const Real dp_values[] = {0.15, 0.1, 0.075};
    const AphiDivFreeValidationFieldKind field_kinds[] = {
        AphiDivFreeValidationFieldKind::Linear2D,
        AphiDivFreeValidationFieldKind::Az2D,
        AphiDivFreeValidationFieldKind::CrossSine3D,
        AphiDivFreeValidationFieldKind::Sinusoidal3DCurlPsi,
    };

    size_t hard_pass_count = 0;
    size_t hard_pass_expected = 0;

    for (const Real dp_0 : dp_values)
    {
        const Real core_shell = core_shell_scale * dp_0;
        for (const AphiDivFreeValidationFieldKind field_kind : field_kinds)
        {
            const AphiDivFreeValidationDivADiagnosticMetrics metrics = runDivFreeValidationDivADiagnostic(
                ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, field_kind);
            printDivFreeValidationDivADiagnosticMetrics("test_3d_aphi_ck_pairwise_diva_3d_diagnostic", metrics);

            if (std::abs(dp_0 - Real(0.1)) > 1.0e-6)
            {
                continue;
            }
            if (field_kind == AphiDivFreeValidationFieldKind::Linear2D)
            {
                hard_pass_expected += 1;
                if (metrics.pairwise.div_a_relative <= max_pairwise_rel_linear2d)
                {
                    hard_pass_count += 1;
                }
            }
            else if (field_kind == AphiDivFreeValidationFieldKind::Az2D)
            {
                hard_pass_expected += 1;
                if (metrics.pairwise.div_a_relative <= max_pairwise_rel_az2d)
                {
                    hard_pass_count += 1;
                }
            }
            else if (field_kind == AphiDivFreeValidationFieldKind::CrossSine3D)
            {
                hard_pass_expected += 1;
                if (metrics.pairwise.div_a_relative <= max_pairwise_rel_crosssine3d)
                {
                    hard_pass_count += 1;
                }
            }
        }
    }

    const AphiDivFreeValidationDivADiagnosticMetrics sinusoidal =
        runDivFreeValidationDivADiagnostic(ac, av, 0.1, body_length, body_height, body_width, 0.25, sigma, nu,
                                           AphiDivFreeValidationFieldKind::Sinusoidal3DCurlPsi);
    const bool sinusoidal_pairwise_high = sinusoidal.pairwise.div_a_relative > 0.20;
    const bool passed = hard_pass_count == hard_pass_expected;

    std::cout << "test_3d_aphi_ck_pairwise_diva_3d_diagnostic hard_pass=" << hard_pass_count << "/"
              << hard_pass_expected << " sinusoidal_pairwise_rel=" << sinusoidal.pairwise.div_a_relative
              << " sinusoidal_B_rel=" << sinusoidal.b_corrected.div_a_relative
              << " sinusoidal_pairwise_high=" << (sinusoidal_pairwise_high ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
