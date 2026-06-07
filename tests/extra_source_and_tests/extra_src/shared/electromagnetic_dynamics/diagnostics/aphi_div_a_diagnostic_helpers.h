#ifndef APHI_DIV_A_DIAGNOSTIC_HELPERS_H
#define APHI_DIV_A_DIAGNOSTIC_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "electromagnetic_dynamics/aphi_pairwise_div_a_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_gradient_divergence_debug_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"
#include "general_gradient.h"
#include "kernel_correction_ck.h"

#include <cmath>
#include <functional>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Post-solve div(A) diagnostic discretization (not used by production GMRES apply). */
enum class AphiDivADiagnosticMode
{
    /** Default: matches operator / PC pairwise divA. */
    PairwiseUncorrected,
    /** Legacy B-corrected trace(gradA) for comparison only. */
    BCorrectedTrace
};

/** Post-solve div(A) diagnostic metrics. Not used by production GMRES. */
struct AphiDivAReductionMetrics
{
    Real div_a_real_L2 = 0.0;
    Real div_a_imag_L2 = 0.0;
    Real div_a_L2 = 0.0;
    Real div_a_Linf = 0.0;
    Real grad_a_L2 = 0.0;
    Real div_a_relative = 0.0;
    size_t particle_count = 0;
};

inline bool isFiniteDivAGradEntry(Real div_value, const Matd &grad_matrix)
{
    return std::isfinite(div_value) && grad_matrix.allFinite();
}

inline void accumulateDivAReductionFromParticles(
    BaseParticles &particles, const AphiVariableNames &names, const Vecd *positions, size_t total_real_particles,
    const std::function<bool(const Vecd &)> &in_region, Real &div_real_squared, Real &div_imag_squared,
    Real &grad_squared, Real &div_linf, size_t &particle_count, AphiDivADiagnosticMode mode)
{
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_real);
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Real *div_a_real = particles.getVariableDataByName<Real>(names.diagnostic.div_a_real);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(names.diagnostic.div_a_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    if (mode == AphiDivADiagnosticMode::BCorrectedTrace)
    {
        syncVariableToHost<Matd>(particles, names.solution.a_real + "Gradient");
        syncVariableToHost<Matd>(particles, names.solution.a_imag + "Gradient");
    }
    else
    {
        syncAphiBlockToHost(particles, names.solution);
    }

    const Matd *grad_a_real =
        mode == AphiDivADiagnosticMode::BCorrectedTrace
            ? particles.getVariableDataByName<Matd>(names.solution.a_real + "Gradient")
            : nullptr;
    const Matd *grad_a_imag =
        mode == AphiDivADiagnosticMode::BCorrectedTrace
            ? particles.getVariableDataByName<Matd>(names.solution.a_imag + "Gradient")
            : nullptr;
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(names.solution.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(names.solution.a_imag);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_region && !in_region(positions[i]))
        {
            continue;
        }
        if (!std::isfinite(div_a_real[i]) || !std::isfinite(div_a_imag[i]))
        {
            continue;
        }
        if (mode == AphiDivADiagnosticMode::BCorrectedTrace &&
            (!isFiniteDivAGradEntry(div_a_real[i], grad_a_real[i]) ||
             !isFiniteDivAGradEntry(div_a_imag[i], grad_a_imag[i])))
        {
            continue;
        }

        particle_count += 1;
        const Real vol_i = vol[i];
        div_real_squared += vol_i * div_a_real[i] * div_a_real[i];
        div_imag_squared += vol_i * div_a_imag[i] * div_a_imag[i];
        if (mode == AphiDivADiagnosticMode::BCorrectedTrace)
        {
            grad_squared += vol_i * (grad_a_real[i].squaredNorm() + grad_a_imag[i].squaredNorm());
        }
        else
        {
            grad_squared += vol_i * (a_real[i].squaredNorm() + a_imag[i].squaredNorm());
        }
        div_linf = std::max(div_linf, std::abs(div_a_real[i]));
        div_linf = std::max(div_linf, std::abs(div_a_imag[i]));
    }
}

inline AphiDivAReductionMetrics finalizeDivAReductionMetrics(Real div_real_squared, Real div_imag_squared,
                                                             Real grad_squared, Real div_linf, size_t particle_count)
{
    AphiDivAReductionMetrics metrics;
    metrics.particle_count = particle_count;
    metrics.div_a_real_L2 = std::sqrt(div_real_squared);
    metrics.div_a_imag_L2 = std::sqrt(div_imag_squared);
    metrics.div_a_L2 = std::sqrt(div_real_squared + div_imag_squared);
    metrics.grad_a_L2 = std::sqrt(grad_squared);
    metrics.div_a_Linf = div_linf;
    if (particle_count == 0 || !std::isfinite(metrics.grad_a_L2))
    {
        metrics.div_a_relative = 0.0;
        return metrics;
    }
    metrics.div_a_relative = metrics.div_a_L2 / (metrics.grad_a_L2 + TinyReal);
    return metrics;
}

inline void execBodyBCorrectedDivADiagnosticPipeline(SPHBody &body, Inner<> &inner, const AphiVariableNames &names)
{
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(inner);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient(inner, names.solution.a_real);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient(inner, names.solution.a_imag);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_real(
        body, names.solution.a_real + "Gradient", names.diagnostic.div_a_real);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_imag(
        body, names.solution.a_imag + "Gradient", names.diagnostic.div_a_imag);

    linear_correction_matrix.exec();
    a_real_gradient.exec();
    a_imag_gradient.exec();
    div_a_real.exec();
    div_a_imag.exec();
}

inline void execBodyPairwiseDivADiagnosticPipeline(SPHBody &body, Inner<> &inner, const AphiVariableNames &names)
{
    (void)body;
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorDivergenceCK<Inner<>>> div_a_real(
        inner, names.solution.a_real, names.diagnostic.div_a_real);
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseVectorDivergenceCK<Inner<>>> div_a_imag(
        inner, names.solution.a_imag, names.diagnostic.div_a_imag);
    div_a_real.exec();
    div_a_imag.exec();
}

inline void execBodyDivADiagnosticPipeline(SPHBody &body, Inner<> &inner, const AphiVariableNames &names,
                                           AphiDivADiagnosticMode mode = AphiDivADiagnosticMode::PairwiseUncorrected)
{
    if (mode == AphiDivADiagnosticMode::BCorrectedTrace)
    {
        execBodyBCorrectedDivADiagnosticPipeline(body, inner, names);
    }
    else
    {
        execBodyPairwiseDivADiagnosticPipeline(body, inner, names);
    }
}

inline AphiDivAReductionMetrics hostBodyDivAReductionMetrics(
    BaseParticles &particles, const AphiVariableNames &names, const Vecd *positions, size_t total_real_particles,
    const std::function<bool(const Vecd &)> &in_region = {},
    AphiDivADiagnosticMode mode = AphiDivADiagnosticMode::PairwiseUncorrected)
{
    Real div_real_squared = 0.0;
    Real div_imag_squared = 0.0;
    Real grad_squared = 0.0;
    Real div_linf = 0.0;
    size_t particle_count = 0;
    accumulateDivAReductionFromParticles(particles, names, positions, total_real_particles, in_region, div_real_squared,
                                       div_imag_squared, grad_squared, div_linf, particle_count, mode);
    return finalizeDivAReductionMetrics(div_real_squared, div_imag_squared, grad_squared, div_linf, particle_count);
}

inline AphiDivAReductionMetrics combineDivAReductionMetrics(const AphiDivAReductionMetrics &lhs,
                                                            const AphiDivAReductionMetrics &rhs)
{
    if (lhs.particle_count == 0)
    {
        return rhs;
    }
    if (rhs.particle_count == 0)
    {
        return lhs;
    }

    AphiDivAReductionMetrics combined;
    combined.particle_count = lhs.particle_count + rhs.particle_count;
    const Real div_real_squared = lhs.div_a_real_L2 * lhs.div_a_real_L2 + rhs.div_a_real_L2 * rhs.div_a_real_L2;
    const Real div_imag_squared = lhs.div_a_imag_L2 * lhs.div_a_imag_L2 + rhs.div_a_imag_L2 * rhs.div_a_imag_L2;
    const Real grad_squared = lhs.grad_a_L2 * lhs.grad_a_L2 + rhs.grad_a_L2 * rhs.grad_a_L2;
    combined.div_a_real_L2 = std::sqrt(div_real_squared);
    combined.div_a_imag_L2 = std::sqrt(div_imag_squared);
    combined.div_a_L2 = std::sqrt(div_real_squared + div_imag_squared);
    combined.grad_a_L2 = std::sqrt(grad_squared);
    combined.div_a_Linf = std::max(lhs.div_a_Linf, rhs.div_a_Linf);
    combined.div_a_relative = combined.div_a_L2 / (combined.grad_a_L2 + TinyReal);
    if (!std::isfinite(combined.div_a_relative))
    {
        combined.div_a_relative = 0.0;
    }
    return combined;
}

inline AphiDivAReductionMetrics hostTwoBodyContactDivAReductionMetrics(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
    const std::function<bool(const Vecd &)> &in_region = {})
{
    execBodyDivADiagnosticPipeline(case_setup.left_body, case_setup.left_inner(), names);
    execBodyDivADiagnosticPipeline(case_setup.right_body, case_setup.right_inner(), names);

    BaseParticles &left_particles = case_setup.left_body.getBaseParticles();
    BaseParticles &right_particles = case_setup.right_body.getBaseParticles();
    syncVariableToHost<Vecd>(left_particles, "Position");
    syncVariableToHost<Vecd>(right_particles, "Position");

    const AphiDivAReductionMetrics left_metrics = hostBodyDivAReductionMetrics(
        left_particles, names, left_particles.getVariableDataByName<Vecd>("Position"),
        left_particles.TotalRealParticles(), in_region);
    const AphiDivAReductionMetrics right_metrics = hostBodyDivAReductionMetrics(
        right_particles, names, right_particles.getVariableDataByName<Vecd>("Position"),
        right_particles.TotalRealParticles(), in_region);
    return combineDivAReductionMetrics(left_metrics, right_metrics);
}

inline AphiDivAReductionMetrics hostTeam7ContactDivAReductionMetrics(
    AphiTeam7ThreeBodyContactCase &case_setup, const AphiVariableNames &names,
    const std::function<bool(const Vecd &)> &in_region = {})
{
    execBodyDivADiagnosticPipeline(case_setup.air_body, case_setup.air_inner(), names);
    execBodyDivADiagnosticPipeline(case_setup.coil_body, case_setup.coil_inner(), names);
    execBodyDivADiagnosticPipeline(case_setup.plate_body, case_setup.plate_inner(), names);

    Real div_real_squared = 0.0;
    Real div_imag_squared = 0.0;
    Real grad_squared = 0.0;
    Real div_linf = 0.0;
    size_t particle_count = 0;
    for (SPHBody *body_ptr : {&case_setup.air_body, &case_setup.coil_body, &case_setup.plate_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        syncVariableToHost<Vecd>(particles, "Position");
        accumulateDivAReductionFromParticles(particles, names, particles.getVariableDataByName<Vecd>("Position"),
                                           particles.TotalRealParticles(), in_region, div_real_squared, div_imag_squared,
                                           grad_squared, div_linf, particle_count,
                                           AphiDivADiagnosticMode::PairwiseUncorrected);
    }
    return finalizeDivAReductionMetrics(div_real_squared, div_imag_squared, grad_squared, div_linf, particle_count);
}

inline const char *divAGaugeDiagnosticLevel(Real div_a_relative)
{
    if (!std::isfinite(div_a_relative))
    {
        return "invalid";
    }
    if (div_a_relative < 1.0e-2)
    {
        return "good";
    }
    if (div_a_relative < 1.0e-1)
    {
        return "warn";
    }
    return "high_risk";
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_DIV_A_DIAGNOSTIC_HELPERS_H
