/**
 * Stage 10.11-C: supplemental analytic ghost layer on device pairwise divA.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_ghost_buffer_diva_diagnostic_helpers.h"

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
    const Real min_sinusoidal_reduction = 0.05;
    const Real max_az2d_baseline_boundary_grad_den = 0.05;

    const AphiGhostBufferDivAMetrics az2d =
        runGhostBufferDivAMetrics(ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu,
                                  AphiDivFreeValidationFieldKind::Az2D);
    const AphiGhostBufferDivAMetrics crosssine =
        runGhostBufferDivAMetrics(ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu,
                                  AphiDivFreeValidationFieldKind::CrossSine3D);
    const AphiGhostBufferDivAMetrics sinusoidal =
        runGhostBufferDivAMetrics(ac, av, dp_0, body_length, body_height, body_width, core_shell, sigma, nu,
                                  AphiDivFreeValidationFieldKind::Sinusoidal3DCurlPsi);

    printGhostBufferDivAMetrics("test_3d_aphi_ck_ghost_buffer_diva_diagnostic", az2d);
    printGhostBufferDivAMetrics("test_3d_aphi_ck_ghost_buffer_diva_diagnostic", crosssine);
    printGhostBufferDivAMetrics("test_3d_aphi_ck_ghost_buffer_diva_diagnostic", sinusoidal);

    const bool az2d_lattice_sufficient = az2d.baseline_boundary_grad_den.div_a_relative <= max_az2d_baseline_boundary_grad_den;
    const bool crosssine_lattice_sufficient =
        crosssine.baseline_boundary_grad_den.div_a_relative <= max_az2d_baseline_boundary_grad_den;
    const bool sinusoidal_ghost_helps = sinusoidal.boundary_grad_den_reduction >= min_sinusoidal_reduction;
    const bool passed = az2d_lattice_sufficient && crosssine_lattice_sufficient && sinusoidal_ghost_helps;

    std::cout << "test_3d_aphi_ck_ghost_buffer_diva_diagnostic az2d_lattice_sufficient=" << (az2d_lattice_sufficient ? 1 : 0)
              << " crosssine_lattice_sufficient=" << (crosssine_lattice_sufficient ? 1 : 0)
              << " sinusoidal_ghost_helps=" << (sinusoidal_ghost_helps ? 1 : 0)
              << " sinusoidal_baseline_boundary_gradDen_rel=" << sinusoidal.baseline_boundary_grad_den.div_a_relative
              << " sinusoidal_ghost_boundary_gradDen_rel=" << sinusoidal.ghost_boundary_grad_den.div_a_relative
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
