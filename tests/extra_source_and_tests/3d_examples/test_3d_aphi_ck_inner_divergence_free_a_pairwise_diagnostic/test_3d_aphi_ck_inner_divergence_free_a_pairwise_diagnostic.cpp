/**
 * Stage 10.9-A: pairwise vs B-corrected divA on analytic divergence-free A fields (no GMRES).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h"

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
    const Real core_shell = 2.5 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real max_pairwise_rel_linear2d = 0.05;
    const Real max_pairwise_rel_sinusoidal3d = 0.20;

    const AphiDivergenceFreeADivADiagnosticMetrics linear2d = runDivergenceFreeADivADiagnostic(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu,
        AphiDivergenceFreeAFieldKind::Linear2D);
    printDivergenceFreeADivADiagnosticMetrics("test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic", linear2d);

    const AphiDivergenceFreeADivADiagnosticMetrics sinusoidal3d = runDivergenceFreeADivADiagnostic(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu,
        AphiDivergenceFreeAFieldKind::Sinusoidal3DCurlPsi);
    printDivergenceFreeADivADiagnosticMetrics("test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic",
                                                sinusoidal3d);

    const bool linear_ok = linear2d.core_particles > 0 &&
                           linear2d.pairwise.div_a_relative <= max_pairwise_rel_linear2d;
    const bool sinusoidal_pairwise_high = sinusoidal3d.pairwise.div_a_relative > max_pairwise_rel_sinusoidal3d;
    const bool passed = linear_ok;

    std::cout << "test_3d_aphi_ck_inner_divergence_free_a_pairwise_diagnostic linear_ok=" << (linear_ok ? 1 : 0)
              << " sinusoidal_pairwise_high=" << (sinusoidal_pairwise_high ? 1 : 0)
              << " sinusoidal_B_rel=" << sinusoidal3d.b_corrected.div_a_relative
              << " sinusoidal_pairwise_rel=" << sinusoidal3d.pairwise.div_a_relative
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
