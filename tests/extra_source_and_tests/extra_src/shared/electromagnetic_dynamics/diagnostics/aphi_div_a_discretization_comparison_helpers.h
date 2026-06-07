#ifndef APHI_DIV_A_DISCRETIZATION_COMPARISON_HELPERS_H
#define APHI_DIV_A_DISCRETIZATION_COMPARISON_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/diagnostics/aphi_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/diagnostics/aphi_a_gauge_diagnostic_helpers.h"
#include "electromagnetic_dynamics/aphi_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <cmath>
#include <functional>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiDivADiscretizationComparisonMetrics
{
    Real div_a_b_l2 = 0.0;
    Real div_a_pairwise_l2 = 0.0;
    Real div_a_b_vs_pairwise_l2_diff = 0.0;
    Real grad_div_a_b_l2 = 0.0;
    Real grad_div_a_pairwise_l2 = 0.0;
    Real grad_div_a_b_vs_pairwise_l2_diff = 0.0;
    Real energy_sign_b = 0.0;
    Real energy_sign_pairwise = 0.0;
    size_t core_particles = 0;
};

inline Real hostCoreVolWeightedScalarL2(
    BaseParticles &particles, const std::string &field_name, const Vecd *positions, size_t total_real_particles,
    Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Real>(particles, field_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *values = particles.getVariableDataByName<Real>(field_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!std::isfinite(values[i]))
        {
            continue;
        }
        squared += vol[i] * values[i] * values[i];
    }
    return std::sqrt(squared);
}

inline Real hostCoreVolWeightedCombinedDivAL2(
    BaseParticles &particles, const std::string &div_a_real_name, const std::string &div_a_imag_name,
    const Vecd *positions, size_t total_real_particles, Real body_length, Real body_height, Real body_width,
    Real core_shell)
{
    syncVariableToHost<Real>(particles, div_a_real_name);
    syncVariableToHost<Real>(particles, div_a_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *div_a_real = particles.getVariableDataByName<Real>(div_a_real_name);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(div_a_imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!std::isfinite(div_a_real[i]) || !std::isfinite(div_a_imag[i]))
        {
            continue;
        }
        squared += vol[i] * (div_a_real[i] * div_a_real[i] + div_a_imag[i] * div_a_imag[i]);
    }
    return std::sqrt(squared);
}

inline Real hostCoreVolWeightedCombinedDivAL2Diff(
    BaseParticles &particles, const std::string &div_a_real_b, const std::string &div_a_imag_b,
    const std::string &div_a_real_p, const std::string &div_a_imag_p, const Vecd *positions,
    size_t total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Real>(particles, div_a_real_b);
    syncVariableToHost<Real>(particles, div_a_imag_b);
    syncVariableToHost<Real>(particles, div_a_real_p);
    syncVariableToHost<Real>(particles, div_a_imag_p);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *div_b_real = particles.getVariableDataByName<Real>(div_a_real_b);
    const Real *div_b_imag = particles.getVariableDataByName<Real>(div_a_imag_b);
    const Real *div_p_real = particles.getVariableDataByName<Real>(div_a_real_p);
    const Real *div_p_imag = particles.getVariableDataByName<Real>(div_a_imag_p);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        const Real diff_real = div_b_real[i] - div_p_real[i];
        const Real diff_imag = div_b_imag[i] - div_p_imag[i];
        if (!std::isfinite(diff_real) || !std::isfinite(diff_imag))
        {
            continue;
        }
        squared += vol[i] * (diff_real * diff_real + diff_imag * diff_imag);
    }
    return std::sqrt(squared);
}

inline Real hostCoreVolWeightedGradDivAL2(
    BaseParticles &particles, const std::string &grad_div_a_real_name, const std::string &grad_div_a_imag_name,
    const Vecd *positions, size_t total_real_particles, Real body_length, Real body_height, Real body_width,
    Real core_shell)
{
    syncVariableToHost<Vecd>(particles, grad_div_a_real_name);
    syncVariableToHost<Vecd>(particles, grad_div_a_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *grad_div_a_real = particles.getVariableDataByName<Vecd>(grad_div_a_real_name);
    const Vecd *grad_div_a_imag = particles.getVariableDataByName<Vecd>(grad_div_a_imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!grad_div_a_real[i].allFinite() || !grad_div_a_imag[i].allFinite())
        {
            continue;
        }
        squared += vol[i] * (grad_div_a_real[i].squaredNorm() + grad_div_a_imag[i].squaredNorm());
    }
    return std::sqrt(squared);
}

inline Real hostCoreVolWeightedGradDivAL2Diff(
    BaseParticles &particles, const std::string &grad_div_a_real_b, const std::string &grad_div_a_imag_b,
    const std::string &grad_div_a_real_p, const std::string &grad_div_a_imag_p, const Vecd *positions,
    size_t total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Vecd>(particles, grad_div_a_real_b);
    syncVariableToHost<Vecd>(particles, grad_div_a_imag_b);
    syncVariableToHost<Vecd>(particles, grad_div_a_real_p);
    syncVariableToHost<Vecd>(particles, grad_div_a_imag_p);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *grad_b_real = particles.getVariableDataByName<Vecd>(grad_div_a_real_b);
    const Vecd *grad_b_imag = particles.getVariableDataByName<Vecd>(grad_div_a_imag_b);
    const Vecd *grad_p_real = particles.getVariableDataByName<Vecd>(grad_div_a_real_p);
    const Vecd *grad_p_imag = particles.getVariableDataByName<Vecd>(grad_div_a_imag_p);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        const Vecd diff_real = grad_b_real[i] - grad_p_real[i];
        const Vecd diff_imag = grad_b_imag[i] - grad_p_imag[i];
        if (!diff_real.allFinite() || !diff_imag.allFinite())
        {
            continue;
        }
        squared += vol[i] * (diff_real.squaredNorm() + diff_imag.squaredNorm());
    }
    return std::sqrt(squared);
}

inline Real hostGradDivEnergySign(
    BaseParticles &particles, const AphiBlockNames &block_a, const std::string &grad_div_a_real_name,
    const std::string &grad_div_a_imag_name, const std::string &div_a_real_name, const std::string &div_a_imag_name,
    size_t total_real_particles, const std::function<bool(const Vecd &)> &in_region)
{
    syncAphiBlockToHost(particles, block_a);
    syncVariableToHost<Vecd>(particles, grad_div_a_real_name);
    syncVariableToHost<Vecd>(particles, grad_div_a_imag_name);
    syncVariableToHost<Real>(particles, div_a_real_name);
    syncVariableToHost<Real>(particles, div_a_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");

    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block_a.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block_a.a_imag);
    const Vecd *grad_div_a_real = particles.getVariableDataByName<Vecd>(grad_div_a_real_name);
    const Vecd *grad_div_a_imag = particles.getVariableDataByName<Vecd>(grad_div_a_imag_name);
    const Real *div_a_real = particles.getVariableDataByName<Real>(div_a_real_name);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(div_a_imag_name);

    Real inner_a_grad_div = 0.0;
    Real div_norm_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_region && !in_region(positions[i]))
        {
            continue;
        }
        if (!std::isfinite(div_a_real[i]) || !std::isfinite(div_a_imag[i]) || !grad_div_a_real[i].allFinite() ||
            !grad_div_a_imag[i].allFinite())
        {
            continue;
        }
        const Real vol_i = vol[i];
        inner_a_grad_div +=
            vol_i * (a_real[i].dot(grad_div_a_real[i]) + a_imag[i].dot(grad_div_a_imag[i]));
        div_norm_squared += vol_i * (div_a_real[i] * div_a_real[i] + div_a_imag[i] * div_a_imag[i]);
    }
    const Real inner_a_minus_grad_div = -inner_a_grad_div;
    return div_norm_squared > TinyReal ? inner_a_minus_grad_div / div_norm_squared : 0.0;
}

inline void execBCorrectedDivAGradDivAPipeline(SPHBody &body, Inner<> &inner, const AphiBlockNames &input_block)
{
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(inner);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient(inner, input_block.a_real);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient(inner, input_block.a_imag);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_real(
        body, input_block.a_real + "Gradient", aphiDivAFieldName(input_block.a_real));
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_imag(
        body, input_block.a_imag + "Gradient", aphiDivAFieldName(input_block.a_imag));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> div_a_real_gradient(
        inner, aphiDivAFieldName(input_block.a_real));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> div_a_imag_gradient(
        inner, aphiDivAFieldName(input_block.a_imag));

    linear_correction_matrix.exec();
    a_real_gradient.exec();
    a_imag_gradient.exec();
    div_a_real.exec();
    div_a_imag.exec();
    div_a_real_gradient.exec();
    div_a_imag_gradient.exec();
}

inline void execPairwiseDivAGradDivAPipeline(SPHBody &body, Inner<> &inner, const AphiBlockNames &input_block)
{
    (void)body;
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorDivergenceCK<Inner<>>> div_a_real(
        inner, input_block.a_real, aphiDivAFieldName(input_block.a_real));
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorDivergenceCK<Inner<>>> div_a_imag(
        inner, input_block.a_imag, aphiDivAFieldName(input_block.a_imag));
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseScalarGradientCK<Inner<>>> div_a_real_gradient(
        inner, aphiDivAFieldName(input_block.a_real), aphiGradDivAFieldName(input_block.a_real));
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseScalarGradientCK<Inner<>>> div_a_imag_gradient(
        inner, aphiDivAFieldName(input_block.a_imag), aphiGradDivAFieldName(input_block.a_imag));

    div_a_real.exec();
    div_a_imag.exec();
    div_a_real_gradient.exec();
    div_a_imag_gradient.exec();
}

inline AphiDivADiscretizationComparisonMetrics hostDivADiscretizationComparisonMetrics(
    BaseParticles &particles, const AphiBlockNames &block_a, const std::string &div_a_real_b,
    const std::string &div_a_imag_b, const std::string &grad_div_a_real_b, const std::string &grad_div_a_imag_b,
    const std::string &div_a_real_p, const std::string &div_a_imag_p, const std::string &grad_div_a_real_p,
    const std::string &grad_div_a_imag_p, size_t total_real_particles, Real body_length, Real body_height,
    Real body_width, Real core_shell)
{
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const auto in_core = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };

    AphiDivADiscretizationComparisonMetrics metrics;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_core(positions[i]))
        {
            metrics.core_particles += 1;
        }
    }

    metrics.div_a_b_l2 = hostCoreVolWeightedCombinedDivAL2(particles, div_a_real_b, div_a_imag_b, positions,
                                                           total_real_particles, body_length, body_height, body_width,
                                                           core_shell);
    metrics.div_a_pairwise_l2 = hostCoreVolWeightedCombinedDivAL2(particles, div_a_real_p, div_a_imag_p, positions,
                                                                  total_real_particles, body_length, body_height,
                                                                  body_width, core_shell);
    metrics.div_a_b_vs_pairwise_l2_diff =
        hostCoreVolWeightedCombinedDivAL2Diff(particles, div_a_real_b, div_a_imag_b, div_a_real_p, div_a_imag_p,
                                              positions, total_real_particles, body_length, body_height, body_width,
                                              core_shell);
    metrics.grad_div_a_b_l2 = hostCoreVolWeightedGradDivAL2(particles, grad_div_a_real_b, grad_div_a_imag_b, positions,
                                                            total_real_particles, body_length, body_height, body_width,
                                                            core_shell);
    metrics.grad_div_a_pairwise_l2 =
        hostCoreVolWeightedGradDivAL2(particles, grad_div_a_real_p, grad_div_a_imag_p, positions, total_real_particles,
                                      body_length, body_height, body_width, core_shell);
    metrics.grad_div_a_b_vs_pairwise_l2_diff =
        hostCoreVolWeightedGradDivAL2Diff(particles, grad_div_a_real_b, grad_div_a_imag_b, grad_div_a_real_p,
                                          grad_div_a_imag_p, positions, total_real_particles, body_length, body_height,
                                          body_width, core_shell);
    metrics.energy_sign_b = hostGradDivEnergySign(particles, block_a, grad_div_a_real_b, grad_div_a_imag_b,
                                                  div_a_real_b, div_a_imag_b, total_real_particles, in_core);
    metrics.energy_sign_pairwise = hostGradDivEnergySign(particles, block_a, grad_div_a_real_p, grad_div_a_imag_p,
                                                         div_a_real_p, div_a_imag_p, total_real_particles, in_core);
    return metrics;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_DIV_A_DISCRETIZATION_COMPARISON_HELPERS_H
