/**
 * Sprint 4.5 Task 1: GMRES workspace contamination via full Contact apply_z_to_w.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_workspace_contamination_helpers.h"

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
    const Real sigma_conductor = 2.0;
    const Real sigma_air = 1.0e-4;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 10.0;
    const Real x_interface = 0.5;
    const Real full_apply_contamination_tol = 1.0e-3;
    const Real block_diagonal_invariance_tol = 1.0e-7;

    AphiVariableNames names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    AphiTwoBodyInterfaceCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    setupTwoBodyInterfaceMmsFields(case_setup, names, sigma_conductor, sigma_air, nu, x_interface);
    RegisterAphiGMRESWorkspaceCK register_left_workspace(case_setup.left_body, 1);
    RegisterAphiGMRESWorkspaceCK register_right_workspace(case_setup.right_body, 1);
    (void)register_left_workspace;
    (void)register_right_workspace;
    case_setup.updateRelations();

    const WorkspaceContaminationMetrics metrics = runWorkspaceContaminationCheck(case_setup, names, options);

    const bool full_apply_shows_contamination = metrics.full_apply_rel_diff > full_apply_contamination_tol;
    const bool block_diagonal_invariant = metrics.block_diagonal_rel_diff < block_diagonal_invariance_tol;
    const bool passed = full_apply_shows_contamination && block_diagonal_invariant &&
                        metrics.reference_output_norm > TinyReal;

    std::cout << "test_3d_aphi_ck_gmres_workspace_contamination"
              << " reference_output_norm=" << metrics.reference_output_norm
              << " full_apply_rel_diff=" << metrics.full_apply_rel_diff
              << " block_diagonal_rel_diff=" << metrics.block_diagonal_rel_diff
              << " full_apply_shows_contamination=" << (full_apply_shows_contamination ? 1 : 0)
              << " block_diagonal_invariant=" << (block_diagonal_invariant ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
