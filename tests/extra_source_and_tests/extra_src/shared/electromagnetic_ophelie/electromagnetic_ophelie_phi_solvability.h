#ifndef ELECTROMAGNETIC_OPHELIE_PHI_SOLVABILITY_H
#define ELECTROMAGNETIC_OPHELIE_PHI_SOLVABILITY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_diagnostics.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline const char *phiRhsOperatorKindName(OpheliePhiRhsOperatorKind kind)
{
    return kind == OpheliePhiRhsOperatorKind::LegacyFlux ? "legacy-flux" : "div-sigma-a";
}

inline Real hostVolWeightedRelativeDifference(BaseParticles &particles, const Real *field_a, const Real *field_b,
                                            size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real diff_l2 = 0.0;
    Real ref_l2 = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real d = field_a[i] - field_b[i];
        diff_l2 += vol[i] * d * d;
        ref_l2 += vol[i] * field_a[i] * field_a[i];
    }
    return std::sqrt(diff_l2) / (std::sqrt(ref_l2) + TinyReal);
}

inline Real hostVolWeightedCosine(BaseParticles &particles, const Real *field_a, const Real *field_b,
                                  size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real dot_ab = 0.0;
    Real norm_a = 0.0;
    Real norm_b = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        dot_ab += vol[i] * field_a[i] * field_b[i];
        norm_a += vol[i] * field_a[i] * field_a[i];
        norm_b += vol[i] * field_b[i] * field_b[i];
    }
    return dot_ab / (std::sqrt(norm_a) * std::sqrt(norm_b) + TinyReal);
}

struct OpheliePhiRhsAlignmentMetrics
{
    Real rhs_div_vs_legacy_vol = 0.0;
    Real rhs_cosine_div_legacy = 0.0;
    /** If legacy ≈ −div: cosine(div, −legacy) → 1 and rel_diff(div, −legacy) → 0. */
    Real rhs_div_vs_neg_legacy_vol = 0.0;
    Real rhs_cosine_div_neg_legacy = 0.0;
};

inline OpheliePhiRhsAlignmentMetrics measurePhiRhsAlignmentMetrics(BaseParticles &particles,
                                                                   const Real *rhs_div, const Real *rhs_legacy,
                                                                   size_t total_real_particles)
{
    OpheliePhiRhsAlignmentMetrics metrics;
    metrics.rhs_div_vs_legacy_vol =
        hostVolWeightedRelativeDifference(particles, rhs_div, rhs_legacy, total_real_particles);
    metrics.rhs_cosine_div_legacy = hostVolWeightedCosine(particles, rhs_div, rhs_legacy, total_real_particles);
    StdVec<Real> neg_legacy(total_real_particles, Real(0));
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        neg_legacy[i] = -rhs_legacy[i];
    }
    metrics.rhs_div_vs_neg_legacy_vol =
        hostVolWeightedRelativeDifference(particles, rhs_div, neg_legacy.data(), total_real_particles);
    metrics.rhs_cosine_div_neg_legacy =
        hostVolWeightedCosine(particles, rhs_div, neg_legacy.data(), total_real_particles);
    return metrics;
}

/** ||RHS_div - RHS_legacy||_vol / ||RHS_div||_vol with same Biot A_src. */
template <class ExecutionPolicy>
inline Real measurePhiRhsDivVsLegacyRelativeDifference(SolidBody &glass_body, Inner<> &inner,
                                                       const OphelieGlassFieldNames &names,
                                                       const OphelieParameters &params)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> rhs_div(n, Real(0));
    StdVec<Real> rhs_legacy(n, Real(0));

    OphelieParameters div_params = params;
    div_params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, div_params);
    hostReadScalarField(particles, names.phi_rhs_imag, rhs_div.data(), n);

    OphelieParameters legacy_params = params;
    legacy_params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, legacy_params);
    hostReadScalarField(particles, names.phi_rhs_imag, rhs_legacy.data(), n);

    return measurePhiRhsAlignmentMetrics(particles, rhs_div.data(), rhs_legacy.data(), n).rhs_div_vs_legacy_vol;
}

/** ||L phi - b||_vol / ||b||_vol for current phi and chosen RHS (LHS already applied to phi_lhs_imag). */
template <class ExecutionPolicy>
inline Real measurePhiEqResVolAtCurrentPhi(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                           const OphelieParameters &params, OpheliePhiRhsOperatorKind rhs_kind)
{
    OphelieParameters rhs_params = params;
    rhs_params.phi_rhs_operator_kind_ = rhs_kind;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, rhs_params);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    return hostPhiEqResVolFromCurrentLhsRhs(glass_body.getBaseParticles(), names,
                                            glass_body.getBaseParticles().TotalRealParticles());
}

struct OpheliePhiRhsSolveDiagnostic
{
    OpheliePhiRhsOperatorKind rhs_kind = OpheliePhiRhsOperatorKind::DivSigmaA;
    Real phi_eq_res_vol = 0.0;
    Real phi_solver_metric = 0.0;
    Real div_j_l2 = 0.0;
    Real div_j_l2_reduction = 0.0;
    /** Same phi: residual if the other RHS were used in the equation. */
    Real cross_eq_res_other_rhs = 0.0;
};

template <class ExecutionPolicy>
inline OpheliePhiRhsSolveDiagnostic runPhiRhsSolveDiagnostic(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters params,
    OpheliePhiRhsOperatorKind rhs_kind, Real div_j_l2_level0, Real characteristic_length)
{
    params.phi_rhs_operator_kind_ = rhs_kind;
    OpheliePhiRhsSolveDiagnostic result;
    result.rhs_kind = rhs_kind;

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, names.phi_imag);
    zero_phi.exec();

    result.phi_solver_metric = solvePhiImag<ExecutionPolicy>(glass_body, inner, names, params);

    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    result.phi_eq_res_vol =
        hostPhiEqResVolFromCurrentLhsRhs(glass_body.getBaseParticles(), names, glass_body.getBaseParticles().TotalRealParticles());
    const OpheliePhiRhsOperatorKind other_rhs =
        rhs_kind == OpheliePhiRhsOperatorKind::DivSigmaA ? OpheliePhiRhsOperatorKind::LegacyFlux
                                                         : OpheliePhiRhsOperatorKind::DivSigmaA;
    result.cross_eq_res_other_rhs =
        measurePhiEqResVolAtCurrentPhi<ExecutionPolicy>(glass_body, inner, names, params, other_rhs);

    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(inner, names);
    StateDynamics<ExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, names, params);
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_grad_phi.exec();
    compute_ejq_with_phi.exec();

    const OphelieDivJMetrics div_j_phi =
        computeOphelieDivJImag<ExecutionPolicy>(glass_body, inner, names, params, characteristic_length);
    result.div_j_l2 = div_j_phi.div_j_weighted_l2;
    result.div_j_l2_reduction = div_j_l2_level0 / (div_j_phi.div_j_weighted_l2 + TinyReal);
    return result;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_SOLVABILITY_H
