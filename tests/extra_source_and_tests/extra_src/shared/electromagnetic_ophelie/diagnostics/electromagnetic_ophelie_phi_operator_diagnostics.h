#ifndef ELECTROMAGNETIC_OPHELIE_PHI_OPERATOR_DIAGNOSTICS_H
#define ELECTROMAGNETIC_OPHELIE_PHI_OPERATOR_DIAGNOSTICS_H

#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_mms_helpers.h"
#include "electromagnetic_ophelie_phi_solvability.h"
#include "update_body_relation.h"

#include <cmath>
#include <iostream>
#include <random>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OpheliePhiOperatorLinearityMetrics
{
    Real repeat_apply_rel = 0.0;
    Real linearity_add_rel = 0.0;
    Real linearity_scale_rel = 0.0;
};

inline Real hostVolWeightedRelativeDifferenceFromVectors(BaseParticles &particles, const Real *field_a,
                                                         const Real *field_b, size_t total_real_particles)
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

template <class ExecutionPolicy>
inline void applyOpheliePhiLhsToHostVector(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                                           const OphelieParameters &params, const Real *phi_values, Real *lhs_values,
                                           size_t total_real_particles)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    hostAssignScalarField(particles, names.phi_imag, phi_values, total_real_particles);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
    hostReadScalarField(particles, names.phi_lhs_imag, lhs_values, total_real_particles);
}

template <class ExecutionPolicy>
inline OpheliePhiOperatorLinearityMetrics evaluateOpheliePhiLhsLinearityMetrics(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    unsigned int seed = 42U)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> phi_x(n, Real(0));
    StdVec<Real> phi_y(n, Real(0));
    StdVec<Real> lhs_x(n, Real(0));
    StdVec<Real> lhs_y(n, Real(0));
    StdVec<Real> lhs_x_repeat(n, Real(0));
    StdVec<Real> lhs_xy(n, Real(0));
    StdVec<Real> lhs_x_plus_y(n, Real(0));
    StdVec<Real> lhs_scaled_x(n, Real(0));
    StdVec<Real> phi_x_plus_y(n, Real(0));
    StdVec<Real> phi_scaled_x(n, Real(0));

    std::mt19937 rng(seed);
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    for (size_t i = 0; i < n; ++i)
    {
        phi_x[i] = dist(rng);
        phi_y[i] = dist(rng);
    }

    const Real alpha = Real(1.75);
    applyOpheliePhiLhsToHostVector<ExecutionPolicy>(glass_body, inner, names, params, phi_x.data(), lhs_x.data(), n);
    applyOpheliePhiLhsToHostVector<ExecutionPolicy>(glass_body, inner, names, params, phi_x.data(), lhs_x_repeat.data(),
                                                    n);
    applyOpheliePhiLhsToHostVector<ExecutionPolicy>(glass_body, inner, names, params, phi_y.data(), lhs_y.data(), n);

    for (size_t i = 0; i < n; ++i)
    {
        phi_x_plus_y[i] = phi_x[i] + phi_y[i];
        phi_scaled_x[i] = alpha * phi_x[i];
    }
    applyOpheliePhiLhsToHostVector<ExecutionPolicy>(glass_body, inner, names, params, phi_x_plus_y.data(),
                                                    lhs_x_plus_y.data(), n);
    applyOpheliePhiLhsToHostVector<ExecutionPolicy>(glass_body, inner, names, params, phi_scaled_x.data(),
                                                    lhs_scaled_x.data(), n);

    for (size_t i = 0; i < n; ++i)
    {
        lhs_xy[i] = lhs_x[i] + lhs_y[i];
    }

    OpheliePhiOperatorLinearityMetrics metrics;
    metrics.repeat_apply_rel =
        hostVolWeightedRelativeDifferenceFromVectors(particles, lhs_x.data(), lhs_x_repeat.data(), n);
    metrics.linearity_add_rel =
        hostVolWeightedRelativeDifferenceFromVectors(particles, lhs_x_plus_y.data(), lhs_xy.data(), n);
    for (size_t i = 0; i < n; ++i)
    {
        lhs_x_repeat[i] = alpha * lhs_x[i];
    }
    metrics.linearity_scale_rel =
        hostVolWeightedRelativeDifferenceFromVectors(particles, lhs_scaled_x.data(), lhs_x_repeat.data(), n);
    return metrics;
}

struct OpheliePhiDiscreteSelfConsistencyMetrics
{
    Real lhs_minus_rhs_vol_l2 = 0.0;
    Real rhs_vol_l2 = 0.0;
    Real discrete_self_consistency_rel = 0.0;
};

/** Manufactured phi_exact -> A_h=-G_h(phi)/omega -> b_h=-omega D_h(sigma A) -> compare L_h(phi) vs b_h. */
template <class ExecutionPolicy>
inline OpheliePhiDiscreteSelfConsistencyMetrics evaluateOpheliePhiDiscreteSelfConsistency(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const Vecd &center, const Vecd &halfsize)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    assignManufacturedPhiExactCosine<ExecutionPolicy>(glass_body, names, center, halfsize);
    assignManufacturedASrcFromPhiExact<ExecutionPolicy>(glass_body, inner, names, params, center, halfsize,
                                                      OpheliePhiMmsSourceKind::DiscreteGrad);

    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, params);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);

    const OpheliePhiEquationResidualMetrics metrics = computeHostPhiEquationResidual(particles, names, n);
    OpheliePhiDiscreteSelfConsistencyMetrics result;
    result.lhs_minus_rhs_vol_l2 = std::sqrt(std::max(metrics.eq_res_vol_l2 * metrics.eq_res_vol_l2 * metrics.rhs_vol_l2,
                                                     Real(0)));
    syncVariableToHost<Real>(particles, names.phi_rhs_imag);
    syncVariableToHost<Real>(particles, names.phi_lhs_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *rhs = particles.getVariableDataByName<Real>(names.phi_rhs_imag);
    const Real *lhs = particles.getVariableDataByName<Real>(names.phi_lhs_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real residual_l2 = 0.0;
    Real rhs_l2 = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real residual = lhs[i] - rhs[i];
        residual_l2 += vol[i] * residual * residual;
        rhs_l2 += vol[i] * rhs[i] * rhs[i];
    }
    result.lhs_minus_rhs_vol_l2 = std::sqrt(residual_l2);
    result.rhs_vol_l2 = std::sqrt(rhs_l2);
    result.discrete_self_consistency_rel = result.lhs_minus_rhs_vol_l2 / (result.rhs_vol_l2 + TinyReal);
    return result;
}

struct OpheliePhiCompatibleOperatorMetrics
{
    Real rhs_unc_vs_compatible_vol = 0.0;
    Real rhs_cosine_unc_compatible = 0.0;
    Real eq_res_vol_uncorrected = 0.0;
    Real eq_res_vol_compatible = 0.0;
};

/** Compare D(sigma A) vs D_c(sigma A) and L_h vs b_h at phi=0 on current A_src. */
template <class ExecutionPolicy>
inline OpheliePhiCompatibleOperatorMetrics evaluateOpheliePhiCompatibleVsUncorrectedOperators(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, OphelieParameters params)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> rhs_unc(n, Real(0));
    StdVec<Real> rhs_comp(n, Real(0));

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, names.phi_imag);
    zero_phi.exec();

    OphelieParameters unc_params = params;
    unc_params.phi_compatible_correction_ = false;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, unc_params);
    hostReadScalarField(particles, names.phi_rhs_imag, rhs_unc.data(), n);

    OphelieParameters comp_params = params;
    comp_params.phi_compatible_correction_ = true;
    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, comp_params);
    hostReadScalarField(particles, names.phi_rhs_imag, rhs_comp.data(), n);

    OpheliePhiCompatibleOperatorMetrics metrics;
    metrics.rhs_unc_vs_compatible_vol =
        hostVolWeightedRelativeDifferenceFromVectors(particles, rhs_unc.data(), rhs_comp.data(), n);
    metrics.rhs_cosine_unc_compatible = hostVolWeightedCosine(particles, rhs_unc.data(), rhs_comp.data(), n);

    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, unc_params);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, unc_params);
    metrics.eq_res_vol_uncorrected = hostPhiEqResVolFromCurrentLhsRhs(particles, names, n);

    setupOpheliePhiImagRhsFromASrc<ExecutionPolicy>(glass_body, inner, names, comp_params);
    applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, comp_params);
    metrics.eq_res_vol_compatible = hostPhiEqResVolFromCurrentLhsRhs(particles, names, n);
    return metrics;
}

inline void logOpheliePhiCompatibleOperatorMetrics(const OpheliePhiCompatibleOperatorMetrics &metrics)
{
    std::cout << "[ophelie] phi_compatible_ops: rhs_unc_vs_compatible_vol=" << metrics.rhs_unc_vs_compatible_vol
              << " rhs_cosine_unc_compatible=" << metrics.rhs_cosine_unc_compatible
              << " eq_res_vol_uncorrected=" << metrics.eq_res_vol_uncorrected
              << " eq_res_vol_compatible=" << metrics.eq_res_vol_compatible << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_OPERATOR_DIAGNOSTICS_H
