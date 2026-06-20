#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_LITERATURE_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_LITERATURE_H

#include "electromagnetic_ophelie_cli.h"
#include "electromagnetic_ophelie_diagnostics.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_gmres.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "electromagnetic_ophelie_phi_operator_diagnostics.h"
#include "electromagnetic_ophelie_phi_boundary_diagnostics.h"
#include "electromagnetic_ophelie_phi_rhs_diagnostics.h"
#include "electromagnetic_ophelie_phi_solvability.h"
#include "electromagnetic_ophelie_edge_flux_diagnostics.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_phi_component.h"
#include "electromagnetic_ophelie_vector_divergence_diagnostics.h"
#include "electromagnetic_ophelie_postprocess.h"
#include "electromagnetic_ophelie_register_fields.h"

#include <cmath>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Jacoutot 2008 OPHELIE-aligned run profile (geometry may still be reduced). */
struct OphelieFrenchLiteratureProfile
{
    bool enabled = false;
    /** Scale coil current so P_raw matches target (Q ∝ I²); no post-hoc field scaling. */
    bool calibrate_coil_current = true;
    /** Minimum acceptable divJ_L2_L0 / divJ_L2_phi (literature target ≥ 1.5). */
    Real div_j_l2_reduction_min = 1.25;
    /** Relative tolerance on P_raw vs target after coil-current calibration. */
    Real power_raw_tolerance = 0.05;
};

inline void applyFrenchLiteratureMode(OphelieParameters &params, OphelieTestCliOptions &cli_options,
                                      OphelieFrenchLiteratureProfile &profile)
{
    profile.enabled = true;
    params.enable_phi_correction_ = true;
    params.enable_power_scaling_ = false;
    cli_options.no_power_scaling = true;
    params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
    params.phi_gmres_max_outer_iterations_ = 120;
    params.phi_gmres_restart_dimension_ = 80;
    params.phi_gmres_eq_res_tolerance_ = 0.0;
    params.phi_eq_res_vol_gate_ = 0.65;
    params.phi_jacobi_max_iterations_ = 8000;
    if (!cli_options.phi_projection_operator_user_set && !cli_options.phi_lhs_operator_user_set &&
        !cli_options.phi_rhs_operator_user_set && !cli_options.ophelie_current_form_user_set)
    {
        applyOpheliePhiProjectionOperatorKind(params, OpheliePhiProjectionOperatorKind::DivGrad);
    }
    params.phi_post_jacobi_refinement_iterations_ = 0;
    std::cout << "[ophelie] literature-mode: phi=ON power_scaling=OFF phi_solver="
              << phiSolverKindName(params.phi_solver_kind_) << " calibrate_current="
              << (profile.calibrate_coil_current ? 1 : 0) << " divJ_L2_red_min=" << profile.div_j_l2_reduction_min
              << std::endl;
}

inline Real calibrateFrenchCoilCurrentToTargetPower(OphelieFrenchReducedCaseParams &french, OphelieParameters &params,
                                                    Real joule_power_raw)
{
    if (params.target_joule_power_ <= TinyReal || joule_power_raw <= TinyReal)
    {
        return 1.0;
    }
    const Real current_scale = std::sqrt(params.target_joule_power_ / joule_power_raw);
    french.coil.current_per_loop *= current_scale;
    french.ampere_turns = french.coil.current_per_loop * static_cast<Real>(french.coil.num_loops);
    syncFrenchReducedToParameters(french, params);
    std::cout << "[ophelie] literature calibrate: P_raw=" << joule_power_raw << " W -> I_per_loop="
              << french.coil.current_per_loop << " A (scale=" << current_scale << ") target_P="
              << params.target_joule_power_ << " W" << std::endl;
    return current_scale;
}

inline void applyOphelieCoilCurrentScale(OphelieFrenchReducedCaseParams &french, OphelieParameters &params)
{
    if (std::abs(params.coil_current_scale_ - 1.0) <= TinyReal)
    {
        return;
    }
    french.coil.current_per_loop *= params.coil_current_scale_;
    french.ampere_turns = french.coil.current_per_loop * static_cast<Real>(french.coil.num_loops);
    syncFrenchReducedToParameters(french, params);
    std::cout << "[ophelie] coil_current_scale=" << params.coil_current_scale_ << " I_per_loop="
              << french.coil.current_per_loop << " ampere_turns=" << french.ampere_turns << std::endl;
}

struct OphelieFrenchEmSolveResult
{
    OphelieDivJMetrics div_j_level0;
    OphelieDivJMetrics div_j_phi;
    Real phi_solver_rel_residual = 0.0;
    Real phi_eq_res_vol = 0.0;
    Real joule_power_raw = 0.0;
    Real joule_power_particle = 0.0;
    /** Graph edge energy diagnostic (not physical Joule power). */
    Real joule_power_graph_edge = 0.0;
    /** Physical Joule power from edge-reconstructed E/J. */
    Real joule_power_recon_edge = 0.0;
    OpheliePhiP0DiagnosticResult phi_p0;
    OphelieEdgeFluxDiagnosticReport edge_flux_report;
    OphelieEdgeFluxPowerMetrics edge_power;
    OphelieEdgeFluxQAntisymMetrics q_antisym;
    OphelieEdgeFluxQSpatialMetrics q_spatial;
    bool edge_flux_report_valid = false;
    bool q_antisym_valid = false;
    bool q_spatial_valid = false;
};

template <class ExecutionPolicy>
inline OphelieFrenchEmSolveResult runFrenchReducedEmPipeline(SolidBody &glass_body, Inner<> &glass_inner,
                                                             const OphelieGlassFieldNames &glass_names,
                                                             OphelieParameters &params,
                                                             const OphelieFrenchReducedCaseParams &french)
{
    OphelieFrenchEmSolveResult result;
    BaseParticles &glass_particles = glass_body.getBaseParticles();
    const size_t n_glass = glass_particles.TotalRealParticles();

    applyMultiloopFilamentBiotToGlass(glass_particles, glass_names, french.coil, params.mu0_, params.softening_length_);

    StateDynamics<ExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_vector_potential(
        glass_body, glass_names);
    combine_vector_potential.exec();

    if (params.phi_biot_divergence_diagnostics_)
    {
        const OphelieBiotSigmaADivergenceDiagnostics biot_div =
            evaluateOphelieBiotSigmaADivergenceDiagnostics<ExecutionPolicy>(
                glass_body, glass_inner, glass_names, french.glass_center, french.glass_radius,
                french.glass_half_height, french.dp);
        logOphelieBiotSigmaADivergenceDiagnostics(biot_div);
        if (!params.phi_biot_divergence_csv_path_.empty())
        {
            appendOphelieBiotSigmaADivergenceCsv(params.phi_biot_divergence_csv_path_, french.dp, biot_div);
        }
    }

    StateDynamics<ExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> compute_ejq_no_phi(glass_body, glass_names,
                                                                                          params);

    if (params.enable_phi_correction_)
    {
        compute_ejq_no_phi.exec();
        result.div_j_level0 =
            computeOphelieDivJImag<ExecutionPolicy>(glass_body, glass_inner, glass_names, params, french.glass_radius);

        const Real boundary_width = params.phi_boundary_distance_factor_ * french.dp;
        if (params.phi_boundary_diagnostics_)
        {
            result.phi_p0.boundary_jn_level0 = computeFrenchCylinderBoundaryJnMetrics(
                glass_particles, glass_names, n_glass, french, params, boundary_width);
        }

        OpheliePhiBoundaryGeometryContext boundary_geom;
        boundary_geom.normal_source = params.phi_boundary_normal_source_;
        boundary_geom.french = french;
        size_t boundary_particle_setup = 0;
        if (params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann)
        {
            boundary_particle_setup =
                setupOpheliePhiBoundaryParticleFields(glass_particles, glass_names, params, boundary_geom, french.dp);
        }

        StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, glass_names.phi_imag);
        zero_phi.exec();
        if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad)
        {
            const OpheliePhiCompatibleOperatorMetrics compatible_ops =
                evaluateOpheliePhiCompatibleVsUncorrectedOperators<ExecutionPolicy>(glass_body, glass_inner,
                                                                                      glass_names, params);
            logOpheliePhiCompatibleOperatorMetrics(compatible_ops);
        }
        setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, glass_inner, glass_names, params);
        OpheliePhiNeumannRhsCorrectionStats neumann_stats;
        finalizeOpheliePhiImagRhsHost(glass_particles, glass_names, params, &boundary_geom, french.dp, &neumann_stats);
        result.phi_p0.neumann_correction_l2 = neumann_stats.correction_l2;
        result.phi_p0.neumann_n_boundary =
            neumann_stats.n_boundary > 0 ? neumann_stats.n_boundary : boundary_particle_setup;
        result.phi_p0.rhs =
            computeOpheliePhiRhsCompatibilityMetrics(glass_particles, glass_names, n_glass);
        result.phi_p0.rhs_zero_mean_projected = params.phi_rhs_project_zero_mean_;
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, glass_inner, glass_names, params);
        result.phi_p0.phi_eq_res_vol_pre_solve =
            hostPhiEqResVolFromCurrentLhsRhs(glass_particles, glass_names, n_glass);
        result.phi_p0.phi_eq_res_vol_pre_solve_after_rhs_proj = result.phi_p0.phi_eq_res_vol_pre_solve;

        OphelieEdgeFluxDiagnosticReport edge_flux_report;
        if (params.phi_edge_flux_diagnostics_)
        {
            edge_flux_report.level0 = evaluateOphelieEdgeFluxResidualBestSignAtCurrentPhi<ExecutionPolicy>(
                glass_body, glass_inner, glass_names, params);
        }

        std::cout << "[ophelie] phi_boundary_mode=" << phiBoundaryModeName(params.phi_boundary_mode_)
                  << " lhs_grad_neumann=" << (params.phi_boundary_lhs_grad_neumann_ ? 1 : 0)
                  << " grad_correction=" << (params.phi_gradient_correction_ ? 1 : 0)
                  << " compatible_correction=" << (params.phi_compatible_correction_ ? 1 : 0)
                  << " phi_pre_solve_eq_res_vol=" << result.phi_p0.phi_eq_res_vol_pre_solve;
        if (params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann)
        {
            std::cout << " neumann_rhs_corr_l2=" << neumann_stats.correction_l2
                      << " neumann_n_boundary=" << neumann_stats.n_boundary;
        }
        std::cout << std::endl;
        if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad)
        {
            BaseParticles &glass_particles = glass_body.getBaseParticles();
            const size_t n_align = glass_particles.TotalRealParticles();
            StdVec<Real> rhs_div(n_align, Real(0));
            StdVec<Real> rhs_legacy(n_align, Real(0));
            OphelieParameters div_p = params;
            div_p.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
            setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, glass_inner, glass_names, div_p);
            hostReadScalarField(glass_particles, glass_names.phi_rhs_imag, rhs_div.data(), n_align);
            OphelieParameters legacy_p = params;
            legacy_p.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
            setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, glass_inner, glass_names, legacy_p);
            hostReadScalarField(glass_particles, glass_names.phi_rhs_imag, rhs_legacy.data(), n_align);
            const OpheliePhiRhsAlignmentMetrics align =
                measurePhiRhsAlignmentMetrics(glass_particles, rhs_div.data(), rhs_legacy.data(), n_align);
            setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, glass_inner, glass_names, params);
            finalizeOpheliePhiImagRhsHost(glass_particles, glass_names, params, &boundary_geom, french.dp, nullptr);
            std::cout << "[ophelie] phi_rhs_solvability: rhs_div_vs_legacy_vol=" << align.rhs_div_vs_legacy_vol
                      << " rhs_cosine_div_neg_legacy=" << align.rhs_cosine_div_neg_legacy
                      << " rhs_div_vs_neg_legacy_vol=" << align.rhs_div_vs_neg_legacy_vol
                      << " rhs_operator=" << phiRhsOperatorKindName(params.phi_rhs_operator_kind_) << std::endl;
        }

        const OpheliePhiRhsFingerprint rhs_after_finalize =
            computeOpheliePhiRhsFingerprint(glass_particles, glass_names.phi_rhs_imag, n_glass);
        logOpheliePhiRhsFingerprint("finalize", rhs_after_finalize);

        if (params.phi_compatible_correction_)
        {
            UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list_prep(glass_body);
            UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation_prep(glass_inner);
            update_cell_linked_list_prep.exec();
            update_inner_relation_prep.exec();
            execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, glass_inner, glass_names, params,
                                                                      Real(1.0e-6), true);
        }

        result.phi_solver_rel_residual =
            solvePhiImagWithCurrentRhs<ExecutionPolicy>(glass_body, glass_inner, glass_names, params);

        if (ophelieUseEdgeFluxElectromotiveRhs(params) && params.edge_flux_complex_)
        {
            const OphelieEdgeFluxComponent real_component = makeOphelieEdgeFluxRealComponent(glass_names, params);
            StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi_real(glass_body, glass_names.phi_real);
            zero_phi_real.exec();
            setupOphelieEdgeFluxComponentRhsFromA<ExecutionPolicy>(glass_body, glass_inner, glass_names, real_component,
                                                                   params);
            finalizeOpheliePhiComponentRhsHost<ExecutionPolicy>(glass_particles, real_component.phi_rhs_field, params);
            const Real real_solver_rel_residual = solvePhiComponentWithCurrentRhs<ExecutionPolicy>(
                glass_body, glass_inner, glass_names, real_component, params);
            applyOpheliePhiComponentLhsOperator<ExecutionPolicy>(glass_body, glass_inner, glass_names, real_component,
                                                                 params);
            std::cout << "[ophelie] edge_flux_complex real_chain: phi_solver_rel_residual=" << real_solver_rel_residual
                      << std::endl;
        }

        const OpheliePhiRhsFingerprint rhs_after_solve =
            computeOpheliePhiRhsFingerprint(glass_particles, glass_names.phi_rhs_imag, n_glass);
        logOpheliePhiRhsFingerprint("post_solve", rhs_after_solve);
        if (!opheliePhiRhsFingerprintsMatch(rhs_after_finalize, rhs_after_solve))
        {
            std::cout << "[ophelie] WARNING: phi_rhs changed between finalize and post_solve" << std::endl;
        }
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, glass_inner, glass_names, params);
        if (params.phi_edge_flux_diagnostics_)
        {
            edge_flux_report.post_phi = evaluateOphelieEdgeFluxResidualBestSignAtCurrentPhi<ExecutionPolicy>(
                glass_body, glass_inner, glass_names, params);
            edge_flux_report.edge_res_red_l2 =
                edge_flux_report.level0.edge_res_l2 / (edge_flux_report.post_phi.edge_res_l2 + TinyReal);
            logOphelieEdgeFluxDiagnosticReport(edge_flux_report);
            if (!params.phi_edge_flux_csv_path_.empty())
            {
                appendOphelieEdgeFluxDiagnosticCsv(params.phi_edge_flux_csv_path_, opheliePhiP0CaseLabel(params),
                                                   french.dp, params, edge_flux_report);
            }
            result.edge_flux_report = edge_flux_report;
            result.edge_flux_report_valid = true;
        }
        if (ophelieUseEdgeFluxElectromotiveRhs(params))
        {
            result.edge_power =
                execOphelieEdgeFluxPostPhiPipeline<ExecutionPolicy>(glass_body, glass_inner, glass_names, params);
            logOphelieEdgeFluxPowerMetrics(result.edge_power);
            result.joule_power_graph_edge = result.edge_power.p_graph_edge;
            result.joule_power_recon_edge = result.edge_power.p_total_recon;
            result.q_antisym = evaluateOphelieEdgeFluxQAntisymmetryForComponent<ExecutionPolicy>(
                glass_body, glass_inner, glass_names, makeOphelieEdgeFluxImagComponent(glass_names, params), params);
            logOphelieEdgeFluxQAntisymMetrics(result.q_antisym, "imag");
            if (params.edge_flux_complex_)
            {
                const OphelieEdgeFluxQAntisymMetrics real_q_antisym =
                    evaluateOphelieEdgeFluxQAntisymmetryForComponent<ExecutionPolicy>(
                        glass_body, glass_inner, glass_names, makeOphelieEdgeFluxRealComponent(glass_names, params),
                        params);
                logOphelieEdgeFluxQAntisymMetrics(real_q_antisym, "real");
            }
            result.q_antisym_valid = true;
            result.q_spatial = computeHostEdgeFluxQSpatialMetrics(
                glass_particles, glass_names, params, n_glass, french.glass_center, Real(0.70) * french.glass_radius,
                Real(0.35) * french.glass_radius, params.q_spatial_max_over_mean_max_,
                params.q_spatial_outer_over_center_min_);
            logOphelieEdgeFluxQSpatialMetrics(result.q_spatial);
            result.q_spatial_valid = true;
        }
        const OpheliePhiEquationResidualMetrics eq_metrics =
            computeHostPhiEquationResidual(glass_particles, glass_names, n_glass);
        result.phi_eq_res_vol = eq_metrics.eq_res_vol_l2;
        result.phi_p0.phi_eq_res_vol_post_solve = eq_metrics.eq_res_vol_l2;
        const OpheliePhiEqResidualShellStats shell_stats = computeHostPhiEqResidualShellStats(
            glass_particles, glass_names, n_glass, french.glass_center,
            Real(0.85) * french.glass_radius);
        std::cout << "[ophelie] phi_eq_res_shell: interior=" << shell_stats.eq_res_interior_vol
                  << " boundary=" << shell_stats.eq_res_boundary_vol
                  << " interior_vol_frac=" << shell_stats.interior_vol_fraction << std::endl;
        if (params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::DivSigmaA)
        {
            const Real cross_legacy = measurePhiEqResVolAtCurrentPhi<ExecutionPolicy>(
                glass_body, glass_inner, glass_names, params, OpheliePhiRhsOperatorKind::LegacyFlux);
            std::cout << "[ophelie] phi_rhs_solvability: cross_eq_res_legacy_at_div_solution=" << cross_legacy
                      << std::endl;
        }
        if (ophelieUseEdgeFluxElectromotiveRhs(params))
        {
            logOphelieEdgeFluxProductionFieldPolicy(params);
            if (params.output_particle_gradient_diagnostics_)
            {
                StateDynamics<ExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_particle_diag(
                    glass_body, glass_names, params);
                UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list_diag(glass_body);
                UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation_diag(glass_inner);
                update_cell_linked_list_diag.exec();
                update_inner_relation_diag.exec();
                execOphelieScalarPhiGradient<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>,
                                             ComputeOphelieScalarPhiGradientCorrectedCK<Inner<>>,
                                             ComputeOpheliePhiGradLinearCorrectionMatrixCK<Inner<>>>(
                    glass_body, glass_inner, glass_names, params);
                applyOpheliePhiBoundaryGradNeumannProjectionDynamics<ExecutionPolicy>(
                    glass_body, glass_names, params, false, &boundary_geom, french.dp);
                compute_ejq_particle_diag.exec();
                StateDynamics<ExecutionPolicy, CopyOpheliePrimaryEJQToParticleDiagnosticCK> copy_particle_diag(
                    glass_body, glass_names);
                copy_particle_diag.exec();
                result.div_j_phi = computeOphelieDivJImag<ExecutionPolicy>(glass_body, glass_inner, glass_names,
                                                                           params, french.glass_radius);
                result.phi_p0.div_j_l2_reduction =
                    result.div_j_level0.div_j_weighted_l2 / (result.div_j_phi.div_j_weighted_l2 + TinyReal);
                std::cout << "[ophelie] edge_flux particle_divJ_diagnostic: divJ_L2_red="
                          << result.phi_p0.div_j_l2_reduction << std::endl;
            }
            syncOphelieEdgeReconToPrimaryEJQ<ExecutionPolicy>(glass_body, glass_names, params);
        }
        else
        {
            StateDynamics<ExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, glass_names,
                                                                                            params);
            UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
            UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(glass_inner);
            update_cell_linked_list.exec();
            update_inner_relation.exec();
            execOphelieScalarPhiGradient<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>,
                                         ComputeOphelieScalarPhiGradientCorrectedCK<Inner<>>,
                                         ComputeOpheliePhiGradLinearCorrectionMatrixCK<Inner<>>>(
                glass_body, glass_inner, glass_names, params);
            applyOpheliePhiBoundaryGradNeumannProjectionDynamics<ExecutionPolicy>(glass_body, glass_names, params,
                                                                                   false, &boundary_geom, french.dp);
            compute_ejq_with_phi.exec();
            result.div_j_phi =
                computeOphelieDivJImag<ExecutionPolicy>(glass_body, glass_inner, glass_names, params, french.glass_radius);
            result.phi_p0.div_j_l2_reduction =
                result.div_j_level0.div_j_weighted_l2 / (result.div_j_phi.div_j_weighted_l2 + TinyReal);
        }
        if (params.phi_boundary_diagnostics_)
        {
            result.phi_p0.boundary_jn_post_phi = computeFrenchCylinderBoundaryJnMetrics(
                glass_particles, glass_names, n_glass, french, params, boundary_width);
        }
        logOpheliePhiP0Diagnostics(result.phi_p0, params);
    }
    else
    {
        compute_ejq_no_phi.exec();
        result.div_j_level0 =
            computeOphelieDivJImag<ExecutionPolicy>(glass_body, glass_inner, glass_names, params, french.glass_radius);
    }

    result.joule_power_particle = hostVolWeightedSum(glass_particles, glass_names.joule_heat, n_glass);
    if (ophelieUseEdgeFluxElectromotiveRhs(params) && result.joule_power_recon_edge > 0.0)
    {
        result.joule_power_raw = result.joule_power_recon_edge;
    }
    else
    {
        result.joule_power_raw = result.joule_power_particle;
    }
    result.phi_p0.joule_power_raw = result.joule_power_raw;
    return result;
}

struct OphelieEdgeFluxLiteratureAcceptance
{
    bool edge_res_red_ok = false;
    bool power_recon_finite = false;
    bool graph_power_diagnostic_ok = false;
    bool graph_power_matches_recon_warning = false;
    bool phi_residual_ok = false;
    bool power_raw_ok = false;
    bool fields_ok = false;
    bool q_antisym_ok = false;
    bool q_spatial_ok = false;
    bool edge_stage1_residual_passed = false;
    bool production_literature_passed = false;
    Real edge_res_red_l2 = 0.0;
    Real p_graph_edge = 0.0;
    Real p_total_recon = 0.0;
    Real p_graph_over_recon = 0.0;
    Real q_antisym_rel_l2 = 0.0;
};

inline OphelieEdgeFluxLiteratureAcceptance evaluateEdgeFluxLiteratureAcceptance(
    const OphelieParameters &params, const OphelieFrenchLiteratureProfile &profile, const OphelieRunMetrics &metrics,
    const OphelieEdgeFluxDiagnosticReport &edge_report, const OphelieEdgeFluxPowerMetrics &power_metrics,
    Real phi_solver_rel_residual, Real joule_power_recon, Real joule_power_particle,
    const OphelieEdgeFluxQAntisymMetrics &q_antisym, bool q_antisym_valid,
    const OphelieEdgeFluxQSpatialMetrics &q_spatial, bool q_spatial_valid)
{
    (void)joule_power_particle;
    OphelieEdgeFluxLiteratureAcceptance acceptance;
    acceptance.edge_res_red_l2 = edge_report.edge_res_red_l2;
    acceptance.p_graph_edge = power_metrics.p_graph_edge;
    acceptance.p_total_recon = joule_power_recon;
    acceptance.p_graph_over_recon = power_metrics.p_graph_over_recon;
    acceptance.edge_res_red_ok = edge_report.edge_res_red_l2 > params.edge_res_red_min_;
    acceptance.power_recon_finite = std::isfinite(joule_power_recon) && joule_power_recon > 0.0;
    acceptance.graph_power_diagnostic_ok =
        std::isfinite(power_metrics.p_graph_edge) && power_metrics.p_graph_edge > 0.0;
    acceptance.graph_power_matches_recon_warning =
        power_metrics.p_graph_over_recon >= params.edge_power_over_recon_soft_min_ &&
        power_metrics.p_graph_over_recon <= params.edge_power_over_recon_soft_max_;
    acceptance.phi_residual_ok =
        params.enable_phi_correction_ &&
        (phi_solver_rel_residual <
         10.0 * (params.phi_solver_kind_ == OpheliePhiSolverKind::GMRES
                     ? params.phi_gmres_tolerance_
                     : params.phi_solver_kind_ == OpheliePhiSolverKind::PCG ? params.phi_pcg_tolerance_
                                                                             : params.phi_jacobi_tolerance_));
    if (profile.calibrate_coil_current && params.target_joule_power_ > TinyReal)
    {
        const Real rel_err = std::abs(joule_power_recon - params.target_joule_power_) / params.target_joule_power_;
        acceptance.power_raw_ok = rel_err <= profile.power_raw_tolerance;
    }
    else
    {
        acceptance.power_raw_ok = joule_power_recon > 0.0;
    }
    acceptance.fields_ok = metrics.n_glass > 0 && metrics.max_a_src > 0.0 && metrics.max_b_src > 0.0 &&
                           power_metrics.joule_heat_edge_recon_max > 0.0 && metrics.max_joule_heat > 0.0;
    acceptance.q_antisym_rel_l2 = q_antisym.q_antisym_rel_l2;
    acceptance.q_antisym_ok =
        q_antisym_valid && q_antisym.q_nonfinite_count == 0 &&
        q_antisym.q_antisym_rel_l2 < params.q_antisym_rel_l2_max_;
    acceptance.q_spatial_ok = !q_spatial_valid || q_spatial.soft_gate_passed;
    acceptance.edge_stage1_residual_passed = acceptance.edge_res_red_ok && acceptance.power_recon_finite &&
                                           acceptance.phi_residual_ok && acceptance.fields_ok &&
                                           acceptance.q_antisym_ok && acceptance.q_spatial_ok &&
                                           !params.enable_self_induction_;
    acceptance.production_literature_passed =
        acceptance.edge_stage1_residual_passed && acceptance.power_raw_ok;
    return acceptance;
}

inline void logEdgeFluxLiteratureAcceptance(const OphelieEdgeFluxLiteratureAcceptance &acceptance)
{
    std::cout << "[ophelie] edge_flux_literature_acceptance: edge_res_red=" << acceptance.edge_res_red_l2
              << " edge_res_gate=" << (acceptance.edge_res_red_ok ? 1 : 0)
              << " P_graph_edge=" << acceptance.p_graph_edge << " P_total_recon=" << acceptance.p_total_recon
              << " P_graph_over_recon=" << acceptance.p_graph_over_recon
              << " P_recon_finite=" << (acceptance.power_recon_finite ? 1 : 0)
              << " P_graph_diagnostic=" << (acceptance.graph_power_diagnostic_ok ? 1 : 0)
              << " P_graph_over_recon_warning=" << (acceptance.graph_power_matches_recon_warning ? 1 : 0)
              << " phi_res=" << (acceptance.phi_residual_ok ? 1 : 0)
              << " P_raw=" << (acceptance.power_raw_ok ? 1 : 0) << " fields=" << (acceptance.fields_ok ? 1 : 0)
              << " q_antisym_rel_l2=" << acceptance.q_antisym_rel_l2
              << " q_antisym_gate=" << (acceptance.q_antisym_ok ? 1 : 0)
              << " q_spatial_gate=" << (acceptance.q_spatial_ok ? 1 : 0)
              << " edge_stage1_residual_passed=" << (acceptance.edge_stage1_residual_passed ? 1 : 0)
              << " production_literature_passed=" << (acceptance.production_literature_passed ? 1 : 0) << std::endl;
}

struct OphelieFrenchLiteratureAcceptance
{
    bool phi_enabled = false;
    bool power_scaling_off = false;
    bool phi_residual_ok = false;
    bool power_raw_ok = false;
    bool div_j_l2_improved = false;
    bool fields_ok = false;
    bool production_operator = false;
    bool passed = false;
    bool production_literature_passed = false;
};

inline OphelieFrenchLiteratureAcceptance evaluateFrenchLiteratureAcceptance(
    const OphelieParameters &params, const OphelieFrenchLiteratureProfile &profile, const OphelieRunMetrics &metrics,
    Real div_j_l2_reduction, Real joule_power_raw, Real phi_eq_res_vol = 0.0)
{
    OphelieFrenchLiteratureAcceptance acceptance;
    acceptance.phi_enabled = params.enable_phi_correction_;
    acceptance.power_scaling_off = !params.enable_power_scaling_;
    if (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad)
    {
        /** DivSigmaGrad on real Biot A: eq_res_vol gate (tightened as post-Jacobi improves). */
        acceptance.phi_residual_ok = params.enable_phi_correction_ && phi_eq_res_vol < params.phi_eq_res_vol_gate_;
    }
    else
    {
        acceptance.phi_residual_ok =
            params.enable_phi_correction_ &&
            (metrics.phi_solver_rel_residual <
             10.0 * (params.phi_solver_kind_ == OpheliePhiSolverKind::GMRES
                         ? params.phi_gmres_tolerance_
                         : params.phi_solver_kind_ == OpheliePhiSolverKind::PCG ? params.phi_pcg_tolerance_
                                                                                 : params.phi_jacobi_tolerance_));
    }
    if (profile.calibrate_coil_current && params.target_joule_power_ > TinyReal)
    {
        const Real rel_err = std::abs(joule_power_raw - params.target_joule_power_) / params.target_joule_power_;
        acceptance.power_raw_ok = rel_err <= profile.power_raw_tolerance;
    }
    else
    {
        acceptance.power_raw_ok = joule_power_raw > 0.0;
    }
    acceptance.div_j_l2_improved = div_j_l2_reduction >= profile.div_j_l2_reduction_min;
    acceptance.fields_ok = metrics.n_glass > 0 && metrics.max_a_src > 0.0 && metrics.max_b_src > 0.0 &&
                           metrics.max_e_imag > 0.0 && metrics.max_j_imag > 0.0 && metrics.max_joule_heat > 0.0 &&
                           metrics.min_joule_heat >= 0.0;
    acceptance.production_operator = opheliePhiProjectionIsProduction(params);
    acceptance.passed = acceptance.phi_enabled && acceptance.power_scaling_off && acceptance.phi_residual_ok &&
                        acceptance.power_raw_ok && acceptance.div_j_l2_improved && acceptance.fields_ok &&
                        !params.enable_self_induction_;
    acceptance.production_literature_passed = acceptance.passed && acceptance.production_operator;
    return acceptance;
}

inline void logFrenchLiteratureAcceptance(const OphelieParameters &params,
                                          const OphelieFrenchLiteratureProfile &profile,
                                          const OphelieFrenchLiteratureAcceptance &acceptance, Real div_j_l2_reduction)
{
    std::cout << "[ophelie] literature_acceptance: projection_operator="
              << phiProjectionOperatorKindName(inferOpheliePhiProjectionOperatorKind(params))
              << " projection_operator_status=" << opheliePhiProjectionRouteStatusName(params)
              << " phi=" << (acceptance.phi_enabled ? 1 : 0) << " no_scaling=" << (acceptance.power_scaling_off ? 1 : 0)
              << " phi_res=" << (acceptance.phi_residual_ok ? 1 : 0)
              << " P_raw=" << (acceptance.power_raw_ok ? 1 : 0)
              << " divJ_L2_red=" << div_j_l2_reduction << " (gate>=" << profile.div_j_l2_reduction_min
              << ", target>=1.5) divJ_gate=" << (acceptance.div_j_l2_improved ? 1 : 0)
              << " fields=" << (acceptance.fields_ok ? 1 : 0)
              << " literature_passed=" << (acceptance.passed ? 1 : 0)
              << " production_literature_passed=" << (acceptance.production_literature_passed ? 1 : 0) << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_LITERATURE_H
