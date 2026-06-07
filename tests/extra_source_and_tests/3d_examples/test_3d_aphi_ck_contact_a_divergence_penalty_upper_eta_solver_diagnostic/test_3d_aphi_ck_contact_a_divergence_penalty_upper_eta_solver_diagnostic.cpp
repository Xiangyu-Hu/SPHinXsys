/**
 * Stage 10.12-C5: Contact upper-eta (0.2) GMRES false-convergence diagnostic.
 * Compares recursive vs true residual, polish_sweep=0 vs 1, and polish curve.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_upper_eta_solver_diagnostic_helpers.h"

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
    const Real tolerance = 1.0e-4;
    const Real primary_eta = AphiADivergencePenaltyResearchDefaults::primary_eta_a;
    const Real upper_eta = AphiADivergencePenaltyResearchDefaults::optional_eta_a;
    const UnsignedInt max_polish_sweeps = 5;
    const Real max_upper_eta_global_true_rel = 1.0e-4;
    const char *test_name = "test_3d_aphi_ck_contact_a_divergence_penalty_upper_eta_solver_diagnostic";

    AphiLhsTestBody scale_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    const AphiCoreOperatorScaleMetrics scale_metrics = hostCoreOperatorScaleMetrics(
        scale_body, body_length, body_height, body_width, core_shell, sigma, nu, omega, phi_gauge_penalty, ac, av);

    const AphiContactUpperEtaSolverDiagnosticRow primary_row = runContactDivFreeSolverDiagnosticRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, primary_eta, scale_metrics, tolerance, x_interface, interface_band_half_width,
        AphiDivFreeValidationFieldKind::Az2D, 1);
    const AphiContactUpperEtaSolverDiagnosticRow upper_no_polish = runContactDivFreeSolverDiagnosticRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, upper_eta, scale_metrics, tolerance, x_interface, interface_band_half_width,
        AphiDivFreeValidationFieldKind::Az2D, 0);
    const AphiContactUpperEtaSolverDiagnosticRow upper_with_polish = runContactDivFreeSolverDiagnosticRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, upper_eta, scale_metrics, tolerance, x_interface, interface_band_half_width,
        AphiDivFreeValidationFieldKind::Az2D, 1);

    printContactUpperEtaSolverDiagnosticRow(test_name, primary_row);
    printContactUpperEtaSolverDiagnosticRow(test_name, upper_no_polish);
    printContactUpperEtaSolverDiagnosticRow(test_name, upper_with_polish);

    AphiTwoBodyInterfaceCase polish_case(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    const AphiBlockNames exact_block{"ExactAReal", "ExactAImag", "ExactPhiReal", "ExactPhiImag"};
    registerDivergenceFreeExactBlock(polish_case.left_body.getBaseParticles(), exact_block);
    registerDivergenceFreeExactBlock(polish_case.right_body.getBaseParticles(), exact_block);
    setupTwoBodyDivFreeMmsFields(polish_case, names, sigma, nu, AphiDivFreeValidationFieldKind::Az2D);
    polish_case.updateRelations();

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = true;
    options.a_divergence_penalty = lambdaAFromEtaA(upper_eta, scale_metrics);
    options.a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;
    options.contact_a_divergence_penalty_stencil = AphiContactADivergencePenaltyStencilMode::InnerOnly;

    const AphiMatrixFreeSolverOptions solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 100);
    const AphiMatrixFreeSolverOptions polish_solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 15);
    const StdVec<AphiTwoBodyPolishSweepSnapshot> curve = runTwoBodyCoupledPolishSweepCurve(
        polish_case, names, options, solver_options, body_length, body_height, body_width, core_shell, x_interface,
        interface_band_half_width, max_polish_sweeps, &polish_solver_options);

    for (const AphiTwoBodyPolishSweepSnapshot &snapshot : curve)
    {
        printContactUpperEtaPolishCurvePoint(test_name, upper_eta, snapshot.polish_sweeps, snapshot.global_true_rel,
                                             snapshot.left_continuous_error, snapshot.right_continuous_error);
    }

    const bool primary_ok = primary_row.gmres_strict_passed && primary_row.mms_solution_finite;
    const bool upper_no_polish_ok =
        contactUpperEtaNoPolishSolutionOk(upper_no_polish, max_upper_eta_global_true_rel);
    const bool upper_reported = std::isfinite(upper_no_polish.exact_consistency_defect) &&
                                std::isfinite(upper_with_polish.exact_consistency_defect);
    const bool polish_curve_reported = !curve.empty();
    const bool upper_with_polish_degraded =
        !std::isfinite(upper_with_polish.global_true_rel) ||
        upper_with_polish.global_true_rel > upper_no_polish.global_true_rel * 10.0;
    const bool passed = primary_ok && upper_no_polish_ok && upper_reported && polish_curve_reported;

    std::cout << test_name << " primary_ok=" << (primary_ok ? 1 : 0)
              << " upper_no_polish_ok=" << (upper_no_polish_ok ? 1 : 0)
              << " upper_reported=" << (upper_reported ? 1 : 0)
              << " upper_with_polish_degraded=" << (upper_with_polish_degraded ? 1 : 0)
              << " polish_curve_points=" << curve.size() << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
