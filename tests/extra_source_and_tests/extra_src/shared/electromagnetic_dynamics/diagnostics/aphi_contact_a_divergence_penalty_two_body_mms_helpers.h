#ifndef APHI_CONTACT_A_DIVERGENCE_PENALTY_TWO_BODY_MMS_HELPERS_H
#define APHI_CONTACT_A_DIVERGENCE_PENALTY_TWO_BODY_MMS_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_eta_sweep_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

class ZeroContactJouleGradFieldsCK : public LocalDynamics
{
  public:
    ZeroContactJouleGradFieldsCK(SPHBody &sph_body, const AphiJouleHeatingFieldNames &field_names)
        : LocalDynamics(sph_body),
          dv_grad_phi_real_(particles_->template getVariableByName<Vecd>(field_names.grad_phi_real)),
          dv_grad_phi_imag_(particles_->template getVariableByName<Vecd>(field_names.grad_phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : grad_phi_real_(encloser.dv_grad_phi_real_->DelegatedData(ex_policy)),
              grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            grad_phi_real_[index_i] = Vecd::Zero();
            grad_phi_imag_[index_i] = Vecd::Zero();
        }

      protected:
        Vecd *grad_phi_real_;
        Vecd *grad_phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_grad_phi_real_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
};

inline void execTwoBodyContactJoulePostProcess(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                               Real omega, const AphiJouleHeatingFieldNames &joule_fields)
{
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>, Contact<>>> left_grad_phi(
        DynamicsArgs(case_setup.left_inner(), names.solution, joule_fields),
        DynamicsArgs(case_setup.left_contact(), names.solution, joule_fields));
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>, Contact<>>> right_grad_phi(
        DynamicsArgs(case_setup.right_inner(), names.solution, joule_fields),
        DynamicsArgs(case_setup.right_contact(), names.solution, joule_fields));
    StateDynamics<MainExecutionPolicy, ZeroContactJouleGradFieldsCK> zero_left_grad(case_setup.left_body, joule_fields);
    StateDynamics<MainExecutionPolicy, ZeroContactJouleGradFieldsCK> zero_right_grad(case_setup.right_body, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> left_electric_field(
        case_setup.left_body, omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> right_electric_field(
        case_setup.right_body, omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> left_joule(case_setup.left_body, names.material,
                                                                               joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> right_joule(case_setup.right_body, names.material,
                                                                                 joule_fields);

    zero_left_grad.exec();
    zero_right_grad.exec();
    left_grad_phi.exec();
    right_grad_phi.exec();
    left_electric_field.exec();
    right_electric_field.exec();
    left_joule.exec();
    right_joule.exec();
    case_setup.updateRelations();
}

inline void setupTwoBodyDivFreeMmsFields(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                         Real sigma, Real nu, AphiDivFreeValidationFieldKind field_kind)
{
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(case_setup.left_body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(case_setup.right_body, sigma, nu,
                                                                                   names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(
        case_setup.left_body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(case_setup.right_body, sigma, nu,
                                                                                       names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_left_sigma(case_setup.left_body, sigma,
                                                                                        names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_right_sigma(case_setup.right_body, sigma,
                                                                                         names.material);
    StateDynamics<MainExecutionPolicy, AssignDivFreeValidationAphiFieldCK> assign_left_exact(
        case_setup.left_body, names.solution, field_kind);
    StateDynamics<MainExecutionPolicy, AssignDivFreeValidationAphiFieldCK> assign_right_exact(
        case_setup.right_body, names.solution, field_kind);

    initialize_left.exec();
    initialize_right.exec();
    set_left_material.exec();
    set_right_material.exec();
    assign_left_sigma.exec();
    assign_right_sigma.exec();
    assign_left_exact.exec();
    assign_right_exact.exec();
}

struct AphiTwoBodyCoreElectromagneticObservables
{
    Real joule_power = 0.0;
    Real E_L2 = 0.0;
    Real J_L2 = 0.0;
    Real E_combined_L2 = 0.0;
    Real J_combined_L2 = 0.0;
};

inline AphiTwoBodyCoreElectromagneticObservables hostTwoBodyCoreObservablesFromJouleFields(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_fields,
    Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma)
{
    AphiTwoBodyCoreElectromagneticObservables combined;
    Real e_combined_squared = 0.0;
    Real j_combined_squared = 0.0;

    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        syncVariableToHost<Vecd>(particles, "Position");
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        const AphiRegionalElectromagneticObservables observables =
            hostCoreObservablesFromJouleFields(particles, names, joule_fields, positions, total_real_particles,
                                               body_length, body_height, body_width, core_shell, sigma);
        combined.joule_power += observables.Joule_power;
        combined.E_L2 = std::sqrt(combined.E_L2 * combined.E_L2 + observables.E_L2 * observables.E_L2);
        combined.J_L2 = std::sqrt(combined.J_L2 * combined.J_L2 + observables.J_L2 * observables.J_L2);
        e_combined_squared += std::pow(hostCoreCombinedElectricFieldL2(particles, joule_fields, positions,
                                                                       total_real_particles, body_length, body_height,
                                                                       body_width, core_shell),
                                       2);
        j_combined_squared += std::pow(hostCoreCombinedCurrentDensityL2(particles, joule_fields, positions,
                                                                        total_real_particles, body_length, body_height,
                                                                        body_width, core_shell),
                                        2);
    }

    combined.E_combined_L2 = std::sqrt(e_combined_squared);
    combined.J_combined_L2 = std::sqrt(j_combined_squared);
    return combined;
}

inline Real hostTwoBodyCoreBlockMaxAbsDifference(AphiTwoBodyInterfaceCase &case_setup, const AphiBlockNames &approx_block,
                                                 const AphiBlockNames &exact_block, Real body_length, Real body_height,
                                                 Real body_width, Real core_shell)
{
    Real max_diff = 0.0;
    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        syncVariableToHost<Vecd>(particles, "Position");
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        max_diff = std::max(max_diff, coreBlockMaxAbsDifference(particles, approx_block, exact_block, positions,
                                                                total_real_particles, body_length, body_height,
                                                                body_width, core_shell));
    }
    return max_diff;
}

inline Real hostTwoBodyCoreDivARelative(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                        Real body_length, Real body_height, Real body_width, Real core_shell)
{
    execBodyDivADiagnosticPipeline(case_setup.left_body, case_setup.left_inner(), names,
                                   AphiDivADiagnosticMode::BCorrectedTrace);
    execBodyDivADiagnosticPipeline(case_setup.right_body, case_setup.right_inner(), names,
                                   AphiDivADiagnosticMode::BCorrectedTrace);

    Real max_rel = 0.0;
    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        max_rel = std::max(max_rel, hostCoreDivARelative(particles, names, total_real_particles, body_length,
                                                         body_height, body_width, core_shell));
    }
    return max_rel;
}

inline Real hostTwoBodyRegionalDivARelative(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                            Real body_length, Real body_height, Real body_width, Real core_shell,
                                            const std::function<bool(const Vecd &)> &in_region)
{
    execBodyDivADiagnosticPipeline(case_setup.left_body, case_setup.left_inner(), names,
                                   AphiDivADiagnosticMode::BCorrectedTrace);
    execBodyDivADiagnosticPipeline(case_setup.right_body, case_setup.right_inner(), names,
                                   AphiDivADiagnosticMode::BCorrectedTrace);

    const AphiDivAReductionMetrics metrics = hostTwoBodyContactDivAReductionMetrics(case_setup, names, in_region);
    return metrics.div_a_relative;
}

struct AphiContactDivFreeTwoBodyMmsRow
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real eta_a = 0.0;
    Real lambda_a = 0.0;
    bool converged = false;
    Real global_true_rel = 0.0;
    Real max_bodywise_true_rel = 0.0;
    Real interface_band_rel = 0.0;
    Real block_linf_error = 0.0;
    Real joule_error_vs_exact = 0.0;
    Real E_combined_error_vs_exact = 0.0;
    Real J_combined_error_vs_exact = 0.0;
    Real core_div_a_rel = 0.0;
    Real interface_div_a_rel = 0.0;
    Real boundary_div_a_rel = 0.0;
    Real discrete_mms_defect = 0.0;
    Real initial_residual_norm = 0.0;
    Real exact_consistency_defect = 0.0;
    std::string breakdown_code;
    AphiInterfaceMmsRunMetrics solver_metrics{};
};

inline AphiContactDivFreeTwoBodyMmsRow runContactDivFreeTwoBodyMmsRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, Real x_interface, Real interface_band_half_width,
    AphiDivFreeValidationFieldKind field_kind, UnsignedInt restart_dimension = 50, UnsignedInt max_outer_iterations = 100,
    AphiContactADivergencePenaltyStencilMode penalty_stencil =
        AphiContactADivergencePenaltyStencilMode::InnerOnly)
{
    AphiContactDivFreeTwoBodyMmsRow row;
    row.field_kind = field_kind;
    row.eta_a = eta_a;
    row.lambda_a = eta_a > TinyReal ? lambdaAFromEtaA(eta_a, scale_metrics) : 0.0;

    AphiTwoBodyInterfaceCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    const AphiJouleHeatingFieldNames joule_fields;
    const AphiBlockNames exact_block{"ExactAReal", "ExactAImag", "ExactPhiReal", "ExactPhiImag"};

    RegisterAphiJouleHeatingFieldsCK register_left_joule(case_setup.left_body, joule_fields);
    RegisterAphiJouleHeatingFieldsCK register_right_joule(case_setup.right_body, joule_fields);
    registerDivergenceFreeExactBlock(case_setup.left_body.getBaseParticles(), exact_block);
    registerDivergenceFreeExactBlock(case_setup.right_body.getBaseParticles(), exact_block);

    setupTwoBodyDivFreeMmsFields(case_setup, names, sigma, nu, field_kind);
    case_setup.updateRelations();

    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_exact(case_setup.left_body, exact_block,
                                                                         names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_exact(case_setup.right_body, exact_block,
                                                                         names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> backup_left_solution(case_setup.left_body, names.t,
                                                                             names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> backup_right_solution(case_setup.right_body, names.t,
                                                                              names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> restore_left_solution(case_setup.left_body, names.solution,
                                                                            names.t);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> restore_right_solution(case_setup.right_body, names.solution,
                                                                              names.t);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_exact_to_solution(case_setup.left_body, names.solution,
                                                                                    exact_block);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_exact_to_solution(
        case_setup.right_body, names.solution, exact_block);

    (void)register_left_joule;
    (void)register_right_joule;
    copy_left_exact.exec();
    copy_right_exact.exec();

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = eta_a > TinyReal;
    options.a_divergence_penalty = row.lambda_a;
    options.a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;
    options.contact_a_divergence_penalty_stencil = penalty_stencil;

    AphiMatrixFreeSolverOptions solver_options =
        defaultCoupledContactGMRESOptions(tolerance, restart_dimension, max_outer_iterations);

    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left_exact(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right_exact(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, omega, options);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_lhs_to_rhs(case_setup.left_body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_lhs_to_rhs(case_setup.right_body, names.rhs,
                                                                            names.lhs);
    apply_left_exact.exec();
    apply_right_exact.exec();
    copy_left_lhs_to_rhs.exec();
    copy_right_lhs_to_rhs.exec();
    case_setup.updateRelations();
    apply_left_exact.exec();
    apply_right_exact.exec();
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
    row.exact_consistency_defect = std::sqrt(lhs_rhs_squared) / (std::sqrt(rhs_squared) + TinyReal);

    // C5: block-GS polish destabilizes coupled GMRES at upper η_A; GMRES alone is stable (polish_sweeps=0).
    const UnsignedInt polish_sweeps =
        eta_a >= AphiADivergencePenaltyResearchDefaults::optional_eta_a - TinyReal ? 0 : 1;
    row.solver_metrics = runTwoBodyCoupledContactInterfaceMms(
        case_setup, names, options, solver_options, body_length, body_height, body_width, core_shell, x_interface,
        interface_band_half_width, polish_sweeps);

    row.converged = row.solver_metrics.solver_result.converged;
    row.global_true_rel = row.solver_metrics.global_true_rel;
    row.max_bodywise_true_rel = row.solver_metrics.max_bodywise_true_rel;
    row.interface_band_rel = row.solver_metrics.interface_band_rel;
    row.discrete_mms_defect = row.solver_metrics.discrete_mms_defect;
    row.initial_residual_norm = row.solver_metrics.solver_result.initial_residual_norm;
    row.breakdown_code = row.solver_metrics.solver_result.breakdown_code_name;
    row.block_linf_error = hostTwoBodyCoreBlockMaxAbsDifference(case_setup, names.solution, exact_block, body_length,
                                                                body_height, body_width, core_shell);

    row.core_div_a_rel = hostTwoBodyCoreDivARelative(case_setup, names, body_length, body_height, body_width,
                                                     core_shell);
    const auto in_interface_band = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell) &&
               std::abs(position[0] - x_interface) <= interface_band_half_width + TinyReal;
    };
    const auto in_boundary = [&](const Vecd &position) {
        return !isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    row.interface_div_a_rel =
        hostTwoBodyRegionalDivARelative(case_setup, names, body_length, body_height, body_width, core_shell,
                                        in_interface_band);
    row.boundary_div_a_rel = hostTwoBodyRegionalDivARelative(case_setup, names, body_length, body_height, body_width,
                                                             core_shell, in_boundary);

    execTwoBodyContactJoulePostProcess(case_setup, names, omega, joule_fields);
    const AphiTwoBodyCoreElectromagneticObservables solved_observables =
        hostTwoBodyCoreObservablesFromJouleFields(case_setup, names, joule_fields, body_length, body_height, body_width,
                                                  core_shell, sigma);

    backup_left_solution.exec();
    backup_right_solution.exec();
    copy_left_exact_to_solution.exec();
    copy_right_exact_to_solution.exec();
    case_setup.updateRelations();
    execTwoBodyContactJoulePostProcess(case_setup, names, omega, joule_fields);
    const AphiTwoBodyCoreElectromagneticObservables exact_observables =
        hostTwoBodyCoreObservablesFromJouleFields(case_setup, names, joule_fields, body_length, body_height, body_width,
                                                  core_shell, sigma);
    restore_left_solution.exec();
    restore_right_solution.exec();
    case_setup.updateRelations();

    row.joule_error_vs_exact = relativeMetricChange(exact_observables.joule_power, solved_observables.joule_power);
    row.E_combined_error_vs_exact =
        relativeMetricChange(exact_observables.E_combined_L2, solved_observables.E_combined_L2);
    row.J_combined_error_vs_exact =
        relativeMetricChange(exact_observables.J_combined_L2, solved_observables.J_combined_L2);
    return row;
}

inline void printContactDivFreeTwoBodyMmsRow(const char *test_name, const AphiContactDivFreeTwoBodyMmsRow &row)
{
    std::cout << test_name << " field=" << divFreeValidationFieldName(row.field_kind) << " eta_A=" << row.eta_a
              << " lambda_A=" << row.lambda_a << " exact_consistency_defect=" << row.exact_consistency_defect
              << " initial_residual=" << row.initial_residual_norm
              << " breakdown=" << row.breakdown_code << " discrete_mms_defect=" << row.discrete_mms_defect
              << " converged=" << (row.converged ? 1 : 0)
              << " global_true_rel=" << row.global_true_rel << " max_bodywise_true_rel=" << row.max_bodywise_true_rel
              << " interface_band_rel=" << row.interface_band_rel << " A_error=" << row.block_linf_error
              << " E_combined_error=" << row.E_combined_error_vs_exact
              << " J_combined_error=" << row.J_combined_error_vs_exact << " Joule_error=" << row.joule_error_vs_exact
              << " core_divA=" << row.core_div_a_rel << " interface_divA=" << row.interface_div_a_rel
              << " boundary_divA=" << row.boundary_div_a_rel << std::endl;
}

inline bool contactDivFreeTwoBodyMmsResearchRowReported(const AphiContactDivFreeTwoBodyMmsRow &row,
                                                        Real max_exact_consistency_defect)
{
    return std::isfinite(row.exact_consistency_defect) && std::isfinite(row.initial_residual_norm) &&
           row.exact_consistency_defect <= max_exact_consistency_defect;
}

inline bool contactDivFreeTwoBodyMmsRowPassed(const AphiContactDivFreeTwoBodyMmsRow &row, Real tolerance,
                                              Real max_block_linf_error)
{
    if (!row.converged || !std::isfinite(row.global_true_rel) || row.global_true_rel > tolerance)
    {
        return false;
    }
    if (!std::isfinite(row.max_bodywise_true_rel) || row.max_bodywise_true_rel > 10.0 * tolerance)
    {
        return false;
    }
    return row.block_linf_error <= max_block_linf_error;
}

inline bool contactDivFreeTwoBodyMmsRowStrictPassed(const AphiContactDivFreeTwoBodyMmsRow &row, Real tolerance,
                                                    Real max_block_linf_error, Real max_exact_consistency_defect,
                                                    Real max_joule_error_vs_exact, Real max_E_combined_error_vs_exact,
                                                    Real max_J_combined_error_vs_exact)
{
    return contactDivFreeTwoBodyMmsRowPassed(row, tolerance, max_block_linf_error) &&
           std::isfinite(row.exact_consistency_defect) &&
           row.exact_consistency_defect <= max_exact_consistency_defect &&
           row.joule_error_vs_exact <= max_joule_error_vs_exact &&
           row.E_combined_error_vs_exact <= max_E_combined_error_vs_exact &&
           row.J_combined_error_vs_exact <= max_J_combined_error_vs_exact;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_A_DIVERGENCE_PENALTY_TWO_BODY_MMS_HELPERS_H
