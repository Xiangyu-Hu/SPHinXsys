#ifndef ELECTROMAGNETIC_OPHELIE_DIAGNOSTICS_H
#define ELECTROMAGNETIC_OPHELIE_DIAGNOSTICS_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_phi_gradient.h"
#include "interaction_algorithms_ck.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"
#include "update_body_relation.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline Real hostVolWeightedVecdNormSquared(BaseParticles &particles, const std::string &variable_name,
                                          size_t total_real_particles)
{
    syncVariableToHost<Vecd>(particles, variable_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *values = particles.getVariableDataByName<Vecd>(variable_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        sum += vol[i] * values[i].squaredNorm();
    }
    return sum;
}

struct OphelieDivJMetrics
{
    Real div_j_weighted_l2 = 0.0;
    Real j_weighted_l2 = 0.0;
    Real div_j_rel = 0.0;
};

inline OphelieDivJMetrics computeHostDivJMetrics(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                size_t total_real_particles, Real characteristic_length)
{
    syncVariableToHost<Real>(particles, names.div_j_imag);
    const Real *div_j = particles.getVariableDataByName<Real>(names.div_j_imag);
    OphelieDivJMetrics metrics;
    metrics.div_j_weighted_l2 = hostVolWeightedNorm(particles, div_j, total_real_particles);
    metrics.j_weighted_l2 = std::sqrt(std::max(hostVolWeightedVecdNormSquared(particles, names.j_imag, total_real_particles),
                                                Real(0)));
    const Real length_scale = std::max(characteristic_length, TinyReal);
    metrics.div_j_rel = metrics.div_j_weighted_l2 / (metrics.j_weighted_l2 / length_scale + TinyReal);
    return metrics;
}

template <class ExecutionPolicy, typename InnerRelationType>
inline OphelieDivJMetrics computeOphelieDivJImag(RealBody &glass_body, InnerRelationType &glass_inner,
                                                 const OphelieGlassFieldNames &names,
                                                 const OphelieParameters &params, Real characteristic_length)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(glass_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    if (opheliePhiUseCorrectedDivergence(params))
    {
        execOpheliePhiGradCorrectionMatrixPrep<ExecutionPolicy>(glass_body, glass_inner, names, params);
    }
    execOphelieVecdDivergence<ExecutionPolicy>(glass_inner, names, params, names.j_imag, names.div_j_imag);
    return computeHostDivJMetrics(glass_body.getBaseParticles(), names, glass_body.getBaseParticles().TotalRealParticles(),
                                  characteristic_length);
}

template <class ExecutionPolicy, typename InnerRelationType>
inline OphelieDivJMetrics computeOphelieDivJImag(RealBody &glass_body, InnerRelationType &glass_inner,
                                                 const OphelieGlassFieldNames &names, Real characteristic_length)
{
    OphelieParameters default_params;
    return computeOphelieDivJImag<ExecutionPolicy>(glass_body, glass_inner, names, default_params,
                                                   characteristic_length);
}

struct OpheliePowerScalingFactors
{
    Real power_raw = 0.0;
    Real power_scale = 1.0;
    Real field_scale = 1.0;
    Real effective_current_amplitude = 0.0;
};

inline OpheliePowerScalingFactors computeOpheliePowerScalingFactors(const OphelieParameters &params, Real joule_power_raw)
{
    OpheliePowerScalingFactors factors;
    factors.power_raw = joule_power_raw;
    if (!params.enable_power_scaling_ || params.target_joule_power_ <= TinyReal)
    {
        factors.power_scale = 1.0;
        factors.field_scale = 1.0;
        factors.effective_current_amplitude = params.current_amplitude_;
        return factors;
    }
    factors.power_scale = params.target_joule_power_ / (joule_power_raw + TinyReal);
    factors.field_scale = std::sqrt(std::max(factors.power_scale, Real(0)));
    factors.effective_current_amplitude = params.current_amplitude_ * factors.field_scale;
    return factors;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DIAGNOSTICS_H
