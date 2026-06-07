/**
 * Stage 10.10-B (continued): pairwise divA root-cause diagnostic — core/boundary split,
 * gradA vs A-norm denominator, core_shell sweep.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_root_cause_helpers.h"

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
    const Real dp_0 = 0.1;
    const Real core_shell_scales[] = {2.5, 3.0, 3.5};
    const AphiDivFreeValidationFieldKind field_kinds[] = {
        AphiDivFreeValidationFieldKind::Linear2D,
        AphiDivFreeValidationFieldKind::Az2D,
        AphiDivFreeValidationFieldKind::CrossSine3D,
    };

    AphiPairwiseDivARootCauseMetrics az2d_ref{};
    bool az2d_ref_set = false;
    size_t grad_den_ratio_pass = 0;
    size_t grad_den_ratio_expected = 0;

    for (const Real shell_scale : core_shell_scales)
    {
        const Real core_shell = shell_scale * dp_0;
        for (const AphiDivFreeValidationFieldKind field_kind : field_kinds)
        {
            const AphiPairwiseDivARootCauseMetrics metrics = runPairwiseDivARootCauseMetrics(
                ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu, field_kind);
            printPairwiseDivARootCauseMetrics("test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic", metrics);

            if (field_kind == AphiDivFreeValidationFieldKind::Az2D && std::abs(shell_scale - 2.5) < 1.0e-6)
            {
                az2d_ref = metrics;
                az2d_ref_set = true;
            }

            if (field_kind == AphiDivFreeValidationFieldKind::Az2D ||
                field_kind == AphiDivFreeValidationFieldKind::CrossSine3D)
            {
                grad_den_ratio_expected += 1;
                const Real anorm_rel = metrics.pairwise_a_norm.div_a_relative;
                const Real grad_rel = metrics.pairwise_grad_den.div_a_relative;
                if (anorm_rel > 0.05 && grad_rel < anorm_rel * 0.5)
                {
                    grad_den_ratio_pass += 1;
                }
            }
        }
    }

    const bool az2d_grad_den_explains_anorm =
        az2d_ref_set && az2d_ref.pairwise_a_norm.div_a_relative > 0.05 &&
        az2d_ref.pairwise_grad_den.div_a_relative < az2d_ref.pairwise_a_norm.div_a_relative * 0.5;
    const bool boundary_fraction_reported =
        az2d_ref_set && std::isfinite(az2d_ref.boundary_div_energy_fraction) &&
        az2d_ref.boundary_div_energy_fraction >= 0.0 && az2d_ref.boundary_div_energy_fraction <= 1.0;
    const bool grad_den_ratio_ok = grad_den_ratio_pass == grad_den_ratio_expected;
    const bool passed = az2d_grad_den_explains_anorm && boundary_fraction_reported && grad_den_ratio_ok;

    std::cout << "test_3d_aphi_ck_pairwise_diva_root_cause_diagnostic az2d_Anorm_rel="
              << (az2d_ref_set ? az2d_ref.pairwise_a_norm.div_a_relative : -1.0)
              << " az2d_gradDen_rel="
              << (az2d_ref_set ? az2d_ref.pairwise_grad_den.div_a_relative : -1.0)
              << " az2d_B_gradDen_rel="
              << (az2d_ref_set ? az2d_ref.b_corrected_grad_norm.div_a_relative : -1.0)
              << " az2d_boundary_div_energy_frac="
              << (az2d_ref_set ? az2d_ref.boundary_div_energy_fraction : -1.0)
              << " grad_den_ratio_pass=" << grad_den_ratio_pass << "/" << grad_den_ratio_expected
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
