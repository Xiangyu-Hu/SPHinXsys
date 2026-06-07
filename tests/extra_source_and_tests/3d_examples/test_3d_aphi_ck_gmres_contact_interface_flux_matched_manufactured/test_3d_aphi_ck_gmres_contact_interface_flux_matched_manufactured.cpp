/**
 * Stage 10-contact-4: two-body Contact coupled GMRES flux-matched interface MMS vs monolithic Inner reference.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <cmath>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

inline AphiInterfaceMmsRunMetrics runMonolithicInterfaceMms(AphiLhsTestBody &test_body, const AphiVariableNames &names,
                                                          const AphiLhsAssemblyOptions &options,
                                                          const AphiMatrixFreeSolverOptions &solver_options,
                                                          Real body_length, Real body_height, Real body_width,
                                                          Real core_shell, Real x_interface,
                                                          Real interface_band_half_width)
{
    AphiInterfaceMmsRunMetrics metrics;

    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact_reference(test_body.body, names.r_hat, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_exact(test_body.body, test_body.inner(), names.solution,
                                                               names.lhs, names.material, options.omega, options);
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_solution(test_body.body, test_body.inner(), names.solution,
                                                                names.lhs, names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_true_residual(
        test_body.body, names.true_residual, names.rhs, names.lhs);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                      solver_options);

    apply_exact.exec();
    copy_exact_reference.exec();
    copy_lhs_to_rhs.exec();
    zero_solution.exec();
    test_body.updateRelations();

    metrics.solver_result = solver.solve();

    apply_solution.exec();
    compute_true_residual.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    metrics.discrete_mms_defect = coreDiscreteMmsOperatorDefect(test_body, names, options, names.solution, names.r_hat);
    metrics.continuous_field_relative_error =
        coreRelativeBlockError(particles, names.solution, names.r_hat, positions, total_real_particles, body_length,
                               body_height, body_width, core_shell);
    const Real interface_band_true_residual =
        hostInterfaceBandTrueResidualNorm(particles, names.true_residual, positions, total_real_particles, body_length,
                                          body_height, body_width, core_shell, x_interface, interface_band_half_width);
    metrics.interface_band_rel = interface_band_true_residual / (metrics.solver_result.initial_residual_norm + TinyReal);
    return metrics;
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
    const Real tolerance = 1.0e-5;
    const Real max_discrete_mms_defect_factor = 10.0;
    const Real x_interface = 0.5;
    const Real interface_band_half_width = 2.0 * dp_0;

    AphiVariableNames names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    const AphiMatrixFreeSolverOptions mono_solver_options = defaultMatrixFreeGMRESOptions(tolerance, 50, 30);
    const AphiMatrixFreeSolverOptions contact_solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 100);

    // --- Monolithic inner reference (9E-1) ---
    AphiLhsTestBody mono_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(mono_body.sph_system);
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

    const AphiInterfaceMmsRunMetrics mono_metrics =
        runMonolithicInterfaceMms(mono_body, names, options, mono_solver_options, body_length, body_height, body_width,
                                  core_shell, x_interface, interface_band_half_width);

    // --- Two-body Contact coupled GMRES ---
    AphiTwoBodyInterfaceCase contact_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    setupTwoBodyInterfaceMmsFields(contact_case, names, sigma_conductor, sigma_air, nu, x_interface);
    contact_case.updateRelations();

    const AphiInterfaceMmsRunMetrics contact_metrics = runTwoBodyCoupledContactInterfaceMms(
        contact_case, names, options, contact_solver_options, body_length, body_height, body_width, core_shell,
        x_interface, interface_band_half_width);

    const Real defect_ratio = contact_metrics.discrete_mms_defect / (mono_metrics.discrete_mms_defect + TinyReal);
    const Real continuous_ratio =
        contact_metrics.continuous_field_relative_error / (mono_metrics.continuous_field_relative_error + TinyReal);

    const bool mono_passed = gmresConvergencePassed(mono_metrics.solver_result, tolerance) &&
                             mono_metrics.discrete_mms_defect <= max_discrete_mms_defect_factor * tolerance;
    const bool contact_solver_passed =
        gmresConvergencePassed(contact_metrics.solver_result, tolerance) &&
        contact_metrics.global_true_rel < 10.0 * tolerance;
    const bool contact_mms_passed =
        contact_metrics.discrete_mms_defect <= max_discrete_mms_defect_factor * tolerance &&
        contact_metrics.interface_band_rel < 1.0e-4;
    const bool passed = mono_passed && contact_solver_passed && contact_mms_passed;

    std::cout << "test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured"
              << " mono_outer=" << mono_metrics.solver_result.outer_iteration_count
              << " mono_discrete_defect=" << mono_metrics.discrete_mms_defect
              << " mono_continuous_error=" << mono_metrics.continuous_field_relative_error
              << " mono_interface_band_rel=" << mono_metrics.interface_band_rel
              << " contact_outer=" << contact_metrics.solver_result.outer_iteration_count
              << " contact_discrete_defect=" << contact_metrics.discrete_mms_defect
              << " contact_interface_band_rel=" << contact_metrics.interface_band_rel
              << " contact_global_true_rel=" << contact_metrics.global_true_rel
              << " contact_max_bodywise_true_rel=" << contact_metrics.max_bodywise_true_rel
              << " contact_left_true_rel=" << contact_metrics.left_true_rel
              << " contact_right_true_rel=" << contact_metrics.right_true_rel
              << " contact_left_continuous=" << contact_metrics.left_continuous_error
              << " contact_right_continuous=" << contact_metrics.right_continuous_error
              << " contact_continuous_error=" << contact_metrics.continuous_field_relative_error
              << " contact_converged=" << contact_metrics.solver_result.converged
              << " contact_rel=" << contact_metrics.solver_result.final_relative_residual
              << " contact_true_rel=" << contact_metrics.solver_result.final_true_relative_residual
              << " defect_ratio=" << defect_ratio << " continuous_ratio=" << continuous_ratio
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
