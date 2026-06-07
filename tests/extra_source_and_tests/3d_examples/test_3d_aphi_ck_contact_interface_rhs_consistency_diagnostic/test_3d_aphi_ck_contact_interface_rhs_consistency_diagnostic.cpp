/**
 * Stage 10-contact-4 diagnostic: interface / RHS / operator consistency for two-body Contact MMS.
 *
 * Scenarios:
 *  A) monolithic exact self-consistency: u=exact, ||Au - b|| / ||b||
 *  B) split coupled RHS at assembly (u=exact on both, before zero)
 *  C) split pin r_hat on both after zero: does stored b match A r_hat ?
 *  D) split left=0, right=r_hat (GMRES left initial state)
 *  E) monolithic vs split apply(r_hat) on matched core positions
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_interface_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

void printRegionalMismatch(const char *label, const AphiRegionalLhsRhsMismatch &metrics)
{
    std::cout << label << " all_rel=" << metrics.all_rel << " core_rel=" << metrics.core_rel
              << " interface_band_rel=" << metrics.interface_band_rel
              << " core_away_interface_rel=" << metrics.core_away_interface_rel << std::endl;
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
    const Real phi_gauge_penalty = 100.0;
    const Real x_interface = 0.5;
    const Real interface_band_half_width = 2.0 * dp_0;
    const Real max_self_consistency_rel = 5.0e-4;
    const Real max_core_away_mono_split_lhs = 5.0e-4;

    AphiVariableNames names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    // --- Monolithic reference ---
    AphiLhsTestBody mono_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_mono(mono_body.body, sigma_conductor, nu,
                                                                                  names);
    StateDynamics<MainExecutionPolicy, AssignPiecewiseSigmaHalfSpaceCK> assign_mono_material(
        mono_body.body, x_interface, sigma_conductor, sigma_air, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignInterfaceFluxMatchedPhiFieldsCK> assign_mono_exact(
        mono_body.body, x_interface, sigma_conductor, sigma_air, names.solution);
    initialize_mono.exec();
    assign_mono_material.exec();
    assign_mono_exact.exec();
    mono_body.updateRelations();
    assembleMonolithicRhsFromExact(mono_body, names, options);

    const AphiRegionalLhsRhsMismatch mono_exact_self = measureBodyLhsRhsMismatch(
        mono_body.body, mono_body.inner(), nullptr, names, options, body_length, body_height, body_width, core_shell,
        x_interface, interface_band_half_width, dp_0, names.solution);

    // --- Two-body split ---
    AphiTwoBodyInterfaceCase split_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    setupTwoBodyInterfaceMmsFields(split_case, names, sigma_conductor, sigma_air, nu, x_interface);
    split_case.updateRelations();
    assembleTwoBodyCoupledRhsFromExact(split_case, names, options);

    // Fresh path: assemble -> zero -> pin r_hat (no intermediate apply/measure).
    AphiTwoBodyInterfaceCase split_pin_only_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    setupTwoBodyInterfaceMmsFields(split_pin_only_case, names, sigma_conductor, sigma_air, nu, x_interface);
    split_pin_only_case.updateRelations();
    assembleTwoBodyCoupledRhsFromExact(split_pin_only_case, names, options);
    zeroTwoBodySolution(split_pin_only_case, names);
    split_pin_only_case.updateRelations();
    pinTwoBodySolutionFromRhat(split_pin_only_case, names);
    split_pin_only_case.updateRelations();
    const Real split_pin_only_solution_rhat_max_diff =
        twoBodyMaxAbsSolutionRhatDifference(split_pin_only_case, names);
    const AphiRegionalLhsRhsMismatch split_pin_only_left = measureBodyLhsRhsMismatch(
        split_pin_only_case.left_body, split_pin_only_case.left_inner(), &split_pin_only_case.left_contact(), names,
        options, body_length, body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0,
        names.solution);
    const AphiRegionalLhsRhsMismatch split_pin_only_right = measureBodyLhsRhsMismatch(
        split_pin_only_case.right_body, split_pin_only_case.right_inner(), &split_pin_only_case.right_contact(), names,
        options, body_length, body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0,
        names.solution);

    AphiTwoBodyInterfaceCase split_zero_left_only_case(dp_0, body_length, body_height, body_width, boundary_width, 0,
                                                       av);
    setupTwoBodyInterfaceMmsFields(split_zero_left_only_case, names, sigma_conductor, sigma_air, nu, x_interface);
    split_zero_left_only_case.updateRelations();
    assembleTwoBodyCoupledRhsFromExact(split_zero_left_only_case, names, options);
    zeroLeftSolutionOnly(split_zero_left_only_case, names);
    split_zero_left_only_case.updateRelations();
    const AphiRegionalLhsRhsMismatch split_zero_left_only_left = measureBodyLhsRhsMismatch(
        split_zero_left_only_case.left_body, split_zero_left_only_case.left_inner(),
        &split_zero_left_only_case.left_contact(), names, options, body_length, body_height, body_width, core_shell,
        x_interface, interface_band_half_width, dp_0, names.solution);

    const AphiRegionalLhsRhsMismatch split_exact_at_assembly_left = measureBodyLhsRhsMismatch(
        split_case.left_body, split_case.left_inner(), &split_case.left_contact(), names, options, body_length,
        body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0, names.solution);
    const AphiRegionalLhsRhsMismatch split_exact_at_assembly_right = measureBodyLhsRhsMismatch(
        split_case.right_body, split_case.right_inner(), &split_case.right_contact(), names, options, body_length,
        body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0, names.solution);

    zeroTwoBodySolution(split_case, names);
    split_case.updateRelations();
    pinTwoBodySolutionFromRhat(split_case, names);
    split_case.updateRelations();
    const Real split_pin_solution_rhat_max_diff = twoBodyMaxAbsSolutionRhatDifference(split_case, names);

    const AphiRegionalLhsRhsMismatch split_pin_rhat_left = measureBodyLhsRhsMismatch(
        split_case.left_body, split_case.left_inner(), &split_case.left_contact(), names, options, body_length,
        body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0, names.solution);
    const AphiRegionalLhsRhsMismatch split_pin_rhat_right = measureBodyLhsRhsMismatch(
        split_case.right_body, split_case.right_inner(), &split_case.right_contact(), names, options, body_length,
        body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0, names.solution);

    zeroTwoBodySolution(split_case, names);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> pin_right_only(split_case.right_body, names.r_hat,
                                                                       names.solution);
    pin_right_only.exec();
    syncAphiBlockToDevice(split_case.right_body.getBaseParticles(), names.solution);
    split_case.updateRelations();

    const AphiRegionalLhsRhsMismatch split_left_zero_right_rhat = measureBodyLhsRhsMismatch(
        split_case.left_body, split_case.left_inner(), &split_case.left_contact(), names, options, body_length,
        body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0, names.solution);

    pinTwoBodySolutionFromRhat(split_case, names);
    split_case.updateRelations();
    const AphiSplitMonoLhsComparison mono_vs_split_lhs =
        compareSplitContactLhsToMonolithic(mono_body, split_case, names, options, body_length, body_height, body_width,
                                           core_shell, x_interface, interface_band_half_width, dp_0);

    const bool mono_self_ok = mono_exact_self.all_rel < max_self_consistency_rel;
    const bool split_assembly_left_ok = split_exact_at_assembly_left.all_rel < max_self_consistency_rel;
    const bool split_assembly_right_ok = split_exact_at_assembly_right.all_rel < max_self_consistency_rel;
    const bool split_zero_left_only_ok = split_zero_left_only_left.all_rel < 5.0e-3;
    const bool split_pin_both_fails_as_expected =
        split_pin_only_left.all_rel > 0.1 && split_pin_only_right.all_rel > 0.1;

    // Diagnostic test passes when sanity checks hold and the zero-neighbor-body failure mode is reproduced.
    const bool passed = mono_self_ok && split_assembly_left_ok && split_assembly_right_ok && split_zero_left_only_ok &&
                        split_pin_both_fails_as_expected;

    std::cout << "test_3d_aphi_ck_contact_interface_rhs_consistency_diagnostic" << std::endl;
    printRegionalMismatch("mono_exact_self", mono_exact_self);
    printRegionalMismatch("split_exact_assembly_left", split_exact_at_assembly_left);
    printRegionalMismatch("split_exact_assembly_right", split_exact_at_assembly_right);
    printRegionalMismatch("split_pin_only_left", split_pin_only_left);
    printRegionalMismatch("split_pin_only_right", split_pin_only_right);
    printRegionalMismatch("split_zero_left_only_left", split_zero_left_only_left);
    printRegionalMismatch("split_pin_rhat_left", split_pin_rhat_left);
    printRegionalMismatch("split_pin_rhat_right", split_pin_rhat_right);
    printRegionalMismatch("split_left0_right_rhat_left", split_left_zero_right_rhat);
    std::cout << "split_pin_only_solution_rhat_max_diff=" << split_pin_only_solution_rhat_max_diff << std::endl;
    std::cout << "split_pin_solution_rhat_max_diff=" << split_pin_solution_rhat_max_diff << std::endl;
    std::cout << "mono_vs_split_lhs matched=" << mono_vs_split_lhs.matched_particles
              << " missing=" << mono_vs_split_lhs.missing_particles
              << " core_away_max_abs_diff=" << mono_vs_split_lhs.max_abs_diff
              << " interface_band_matched=" << mono_vs_split_lhs.interface_band_matched
              << " interface_band_max_abs_diff=" << mono_vs_split_lhs.interface_band_max_abs_diff << std::endl;
    std::cout << "sanity mono_self=" << (mono_self_ok ? 1 : 0) << " split_assembly=" << (split_assembly_left_ok && split_assembly_right_ok ? 1 : 0)
              << " split_zero_left_only=" << (split_zero_left_only_ok ? 1 : 0)
              << " split_pin_both_fails=" << (split_pin_both_fails_as_expected ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
