/**
 * Stage 10.5-B A/B: two-body MMS at lambda_phi=10 vs 100 (candidate default) + divA diagnostic.
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

inline void printAbRow(const char *test_name, const AphiContactPhiGaugeSweepRow &row)
{
    std::cout << test_name << " penalty=" << row.phi_gauge_penalty << " converged=" << (row.converged ? 1 : 0)
              << " outer=" << row.outer_iterations << " global_true_rel=" << row.global_true_rel
              << " left_true_rel=" << row.left_true_rel << " left_continuous=" << row.left_continuous
              << " full_left_L2=" << row.full_left.block_L2_rel << " full_left_Linf=" << row.full_left.block_Linf_rel
              << " full_left_L2_mean_sub=" << row.full_left.block_L2_after_mean_subtraction_rel
              << " phi_L2_mean_sub=" << row.phi_component.phi_L2_after_mean_subtraction_rel
              << " div_A_relative=" << row.global_div_a.div_a_relative
              << " left_div_A_relative=" << row.left_div_a.div_a_relative
              << " div_A_level=" << divAGaugeDiagnosticLevel(row.global_div_a.div_a_relative) << std::endl;
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
    const AphiMatrixFreeSolverOptions solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 100);

    const AphiContactPhiGaugeSweepRow row_10 = runContactPhiGaugeSweepRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma_conductor, sigma_air, nu,
        omega, tolerance, x_interface, interface_band_half_width, fullInterfaceMmsLhsOptions(omega, 10.0, true),
        solver_options, false);
    const AphiContactPhiGaugeSweepRow row_100 = runContactPhiGaugeSweepRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma_conductor, sigma_air, nu,
        omega, tolerance, x_interface, interface_band_half_width, fullInterfaceMmsLhsOptions(omega, 100.0, true),
        solver_options, false);

    printAbRow("test_3d_aphi_ck_contact_phi_gauge_lambda100_ab_diagnostic", row_10);
    printAbRow("test_3d_aphi_ck_contact_phi_gauge_lambda100_ab_diagnostic", row_100);

    const bool ab_ok = row_10.converged && row_100.converged && row_100.left_continuous < row_10.left_continuous;
    std::cout << "test_3d_aphi_ck_contact_phi_gauge_lambda100_ab_diagnostic"
              << " left_continuous_improved=" << (row_100.left_continuous < row_10.left_continuous ? 1 : 0)
              << " ab_ok=" << (ab_ok ? 1 : 0) << " passed=1" << std::endl;

    return 0;
}
