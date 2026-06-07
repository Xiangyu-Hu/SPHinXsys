#ifndef APHI_CONTACT_UPPER_ETA_SOLVER_DIAGNOSTIC_HELPERS_H
#define APHI_CONTACT_UPPER_ETA_SOLVER_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_contact_a_divergence_penalty_two_body_mms_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactUpperEtaSolverDiagnosticRow
{
    Real eta_a = 0.0;
    Real lambda_a = 0.0;
    UnsignedInt polish_sweeps = 0;
    Real exact_consistency_defect = 0.0;
    AphiMatrixFreeSolverResult solver_result{};
    Real global_true_rel = 0.0;
    Real discrete_mms_defect = 0.0;
    bool gmres_strict_passed = false;
    bool mms_solution_finite = false;
};

inline Real hostTwoBodyExactConsistencyDefect(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                              const AphiLhsAssemblyOptions &options)
{
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_lhs_to_rhs(case_setup.left_body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_lhs_to_rhs(case_setup.right_body, names.rhs,
                                                                            names.lhs);

    apply_left.exec();
    apply_right.exec();
    copy_left_lhs_to_rhs.exec();
    copy_right_lhs_to_rhs.exec();
    case_setup.updateRelations();
    apply_left.exec();
    apply_right.exec();
    case_setup.updateRelations();

    Real lhs_rhs_squared = 0.0;
    Real rhs_squared = 0.0;
    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        syncAphiBlockToHost(particles, names.lhs);
        syncAphiBlockToHost(particles, names.rhs);
        StateDynamics<MainExecutionPolicy, AphiBlockLinearCombinationCK> form_difference(
            *body_ptr, names.residual, Real(1), Real(-1), names.lhs, names.rhs);
        form_difference.exec();
        lhs_rhs_squared += std::pow(hostBlockNorm(particles, names.residual, total_real_particles), 2);
        rhs_squared += std::pow(hostBlockNorm(particles, names.rhs, total_real_particles), 2);
    }
    return std::sqrt(lhs_rhs_squared) / (std::sqrt(rhs_squared) + TinyReal);
}

inline AphiContactUpperEtaSolverDiagnosticRow runContactDivFreeSolverDiagnosticRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, Real x_interface, Real interface_band_half_width,
    AphiDivFreeValidationFieldKind field_kind, UnsignedInt polish_sweeps, Real max_recursive_true_gap = 1.0e-4)
{
    AphiContactUpperEtaSolverDiagnosticRow row;
    row.eta_a = eta_a;
    row.lambda_a = eta_a > TinyReal ? lambdaAFromEtaA(eta_a, scale_metrics) : 0.0;
    row.polish_sweeps = polish_sweeps;

    AphiTwoBodyInterfaceCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    const AphiBlockNames exact_block{"ExactAReal", "ExactAImag", "ExactPhiReal", "ExactPhiImag"};

    registerDivergenceFreeExactBlock(case_setup.left_body.getBaseParticles(), exact_block);
    registerDivergenceFreeExactBlock(case_setup.right_body.getBaseParticles(), exact_block);
    setupTwoBodyDivFreeMmsFields(case_setup, names, sigma, nu, field_kind);
    case_setup.updateRelations();

    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_exact(case_setup.left_body, exact_block,
                                                                         names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_exact(case_setup.right_body, exact_block,
                                                                         names.solution);
    copy_left_exact.exec();
    copy_right_exact.exec();

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = eta_a > TinyReal;
    options.a_divergence_penalty = row.lambda_a;
    options.a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;
    options.contact_a_divergence_penalty_stencil = AphiContactADivergencePenaltyStencilMode::InnerOnly;

    row.exact_consistency_defect = hostTwoBodyExactConsistencyDefect(case_setup, names, options);

    const AphiMatrixFreeSolverOptions solver_options =
        defaultCoupledContactGMRESOptions(tolerance, 50, 100);
    const AphiInterfaceMmsRunMetrics metrics = runTwoBodyCoupledContactInterfaceMms(
        case_setup, names, options, solver_options, body_length, body_height, body_width, core_shell, x_interface,
        interface_band_half_width, polish_sweeps);

    row.solver_result = metrics.solver_result;
    row.global_true_rel = metrics.global_true_rel;
    row.discrete_mms_defect = metrics.discrete_mms_defect;
    row.gmres_strict_passed = gmresConvergencePassed(metrics.solver_result, tolerance, 10.0, max_recursive_true_gap);
    row.mms_solution_finite =
        std::isfinite(row.global_true_rel) && std::isfinite(row.discrete_mms_defect);
    return row;
}

inline void printContactUpperEtaSolverDiagnosticRow(const char *test_name,
                                                    const AphiContactUpperEtaSolverDiagnosticRow &row)
{
    std::cout << test_name << " eta_A=" << row.eta_a << " lambda_A=" << row.lambda_a
              << " polish_sweeps=" << row.polish_sweeps << " exact_consistency_defect=" << row.exact_consistency_defect
              << " gmres_converged=" << (row.solver_result.converged ? 1 : 0)
              << " recursive_rel=" << row.solver_result.final_relative_residual
              << " true_rel=" << row.solver_result.final_true_relative_residual
              << " recursive_true_gap=" << row.solver_result.final_recursive_true_gap
              << " gmres_strict_passed=" << (row.gmres_strict_passed ? 1 : 0)
              << " global_true_rel=" << row.global_true_rel << " discrete_mms_defect=" << row.discrete_mms_defect
              << " mms_finite=" << (row.mms_solution_finite ? 1 : 0) << std::endl;
}

inline bool contactUpperEtaNoPolishSolutionOk(const AphiContactUpperEtaSolverDiagnosticRow &row,
                                              Real max_global_true_rel = 1.0e-4)
{
    return row.polish_sweeps == 0 && std::isfinite(row.global_true_rel) && row.global_true_rel <= max_global_true_rel &&
           row.mms_solution_finite && row.gmres_strict_passed;
}

inline void printContactUpperEtaPolishCurvePoint(const char *test_name, Real eta_a, UnsignedInt polish_sweeps,
                                                 Real global_true_rel, Real left_continuous_error,
                                                 Real right_continuous_error)
{
    std::cout << test_name << " eta_A=" << eta_a << " polish_sweeps=" << polish_sweeps
              << " global_true_rel=" << global_true_rel << " left_continuous_err=" << left_continuous_error
              << " right_continuous_err=" << right_continuous_error << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_UPPER_ETA_SOLVER_DIAGNOSTIC_HELPERS_H
