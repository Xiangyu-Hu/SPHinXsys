/**
 * Informational: two-body coupled GMRES + incremental block-GS polish sweep curve.
 * Left-body max-norm vs manufactured exact is tracked but does not gate pass/fail.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"

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
    const Real sigma_conductor = 2.0;
    const Real sigma_air = 1.0e-4;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real tolerance = 1.0e-5;
    const Real x_interface = 0.5;
    const Real interface_band_half_width = 2.0 * dp_0;
    const UnsignedInt max_polish_sweeps = 3;

    AphiVariableNames names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    const AphiMatrixFreeSolverOptions coupled_solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 100);
    const AphiMatrixFreeSolverOptions polish_solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 15);

    AphiTwoBodyInterfaceCase contact_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    setupTwoBodyInterfaceMmsFields(contact_case, names, sigma_conductor, sigma_air, nu, x_interface);
    contact_case.updateRelations();

    const StdVec<AphiTwoBodyPolishSweepSnapshot> curve = runTwoBodyCoupledPolishSweepCurve(
        contact_case, names, options, coupled_solver_options, body_length, body_height, body_width, core_shell,
        x_interface, interface_band_half_width, max_polish_sweeps, &polish_solver_options);

    std::cout << "test_3d_aphi_ck_contact_two_body_polish_sweep_diagnostic"
              << " max_polish_sweeps=" << max_polish_sweeps;
    for (const AphiTwoBodyPolishSweepSnapshot &snapshot : curve)
    {
        std::cout << " sweep=" << snapshot.polish_sweeps << " global_true_rel=" << snapshot.global_true_rel
                  << " left_true_rel=" << snapshot.left_true_rel << " right_true_rel=" << snapshot.right_true_rel
                  << " left_continuous=" << snapshot.left_continuous_error
                  << " right_continuous=" << snapshot.right_continuous_error;
    }
    std::cout << " passed=1" << std::endl;

    return 0;
}
