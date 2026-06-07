/**
 * Stage 10.5-B: phi gauge penalty sweep + phi-only operator comparison on two-body Contact MMS.
 * Informational diagnostic — always exits 0 when all rows complete.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_phi_gauge_sweep_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

inline void printSweepRow(const char *test_name, const AphiContactPhiGaugeSweepRow &row)
{
    std::cout << test_name << " penalty=" << row.phi_gauge_penalty
              << " use_penalty=" << (row.use_phi_gauge_penalty ? 1 : 0)
              << " phi_only=" << (row.phi_only_operator ? 1 : 0) << " converged=" << (row.converged ? 1 : 0)
              << " outer=" << row.outer_iterations << " global_true_rel=" << row.global_true_rel
              << " left_true_rel=" << row.left_true_rel << " left_continuous=" << row.left_continuous
              << " full_left_L2=" << row.full_left.block_L2_rel << " full_left_Linf=" << row.full_left.block_Linf_rel
              << " full_left_L2_mean_sub=" << row.full_left.block_L2_after_mean_subtraction_rel
              << " phi_L2=" << row.phi_component.phi_L2_rel << " phi_Linf=" << row.phi_component.phi_Linf_rel
              << " phi_real_mean=" << row.phi_component.phi_real_mean_error
              << " phi_L2_mean_sub=" << row.phi_component.phi_L2_after_mean_subtraction_rel
              << " core_interior_Linf=" << row.core_interior.block_Linf_rel
              << " core_interior_phi_L2_mean_sub=" << row.core_interior_phi.phi_L2_after_mean_subtraction_rel
              << " div_A_relative=" << row.global_div_a.div_a_relative
              << " div_A_Linf=" << row.global_div_a.div_a_Linf
              << " left_div_A_relative=" << row.left_div_a.div_a_relative
              << " right_div_A_relative=" << row.right_div_a.div_a_relative
              << " div_A_level=" << divAGaugeDiagnosticLevel(row.global_div_a.div_a_relative)
              << std::endl;
}

} // namespace

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
    const Real tolerance = 1.0e-5;
    const Real x_interface = 0.5;
    const Real interface_band_half_width = 2.0 * dp_0;
    const Real baseline_penalty = 10.0;
    const AphiMatrixFreeSolverOptions solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 100);

    const Real penalty_values[] = {0.0, 1.0, 10.0, 100.0};
    UnsignedInt completed_rows = 0;

    std::cout << "test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic begin" << std::endl;

    for (const Real penalty : penalty_values)
    {
        const bool use_penalty = penalty > TinyReal;
        const AphiLhsAssemblyOptions options = fullInterfaceMmsLhsOptions(omega, penalty, use_penalty);
        const AphiContactPhiGaugeSweepRow row = runContactPhiGaugeSweepRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma_conductor, sigma_air,
            nu, omega, tolerance, x_interface, interface_band_half_width, options, solver_options, false);
        printSweepRow("test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic", row);
        completed_rows += 1;
    }

    const AphiLhsAssemblyOptions phi_only_baseline =
        phiOnlyLaplacePenaltyLhsOptions(omega, baseline_penalty, true);
    const AphiContactPhiGaugeSweepRow phi_only_row = runContactPhiGaugeSweepRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma_conductor, sigma_air, nu,
        omega, tolerance, x_interface, interface_band_half_width, phi_only_baseline, solver_options, true);
    printSweepRow("test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic", phi_only_row);
    completed_rows += 1;

    const AphiLhsAssemblyOptions phi_only_no_penalty = phiOnlyLaplacePenaltyLhsOptions(omega, 0.0, false);
    const AphiContactPhiGaugeSweepRow phi_only_off_row = runContactPhiGaugeSweepRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma_conductor, sigma_air, nu,
        omega, tolerance, x_interface, interface_band_half_width, phi_only_no_penalty, solver_options, true);
    printSweepRow("test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic", phi_only_off_row);
    completed_rows += 1;

    std::cout << "test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic"
              << " completed_rows=" << completed_rows << " passed=1" << std::endl;

    return 0;
}
