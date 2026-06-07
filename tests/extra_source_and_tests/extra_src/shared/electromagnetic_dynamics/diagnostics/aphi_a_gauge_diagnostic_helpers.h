#ifndef APHI_A_GAUGE_DIAGNOSTIC_HELPERS_H
#define APHI_A_GAUGE_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
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

inline Real hostAVectorVolDotProduct(BaseParticles &particles, const AphiBlockNames &block_a,
                                     const AphiBlockNames &block_b, size_t total_real_particles,
                                     const std::function<bool(const Vecd &)> &in_region = {})
{
    syncAphiBlockToHost(particles, block_a);
    syncAphiBlockToHost(particles, block_b);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");

    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block_a.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block_a.a_imag);
    const Vecd *b_real = particles.getVariableDataByName<Vecd>(block_b.a_real);
    const Vecd *b_imag = particles.getVariableDataByName<Vecd>(block_b.a_imag);

    Real dot = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_region && !in_region(positions[i]))
        {
            continue;
        }
        dot += vol[i] * (a_real[i].dot(b_real[i]) + a_imag[i].dot(b_imag[i]));
    }
    return dot;
}

inline Real hostDivANormSquared(BaseParticles &particles, const AphiVariableNames &names, size_t total_real_particles,
                                const std::function<bool(const Vecd &)> &in_region = {})
{
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_real);
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");

    const Real *div_a_real = particles.getVariableDataByName<Real>(names.diagnostic.div_a_real);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(names.diagnostic.div_a_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    Real squared = 0.0;
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
        squared += vol[i] * (div_a_real[i] * div_a_real[i] + div_a_imag[i] * div_a_imag[i]);
    }
    return squared;
}

struct AphiGradDivSignEnergyMetrics
{
    Real inner_a_plus_graddiv_a = 0.0;
    Real inner_a_minus_graddiv_a = 0.0;
    Real div_a_norm_squared = 0.0;
    Real ratio_plus = 0.0;
    Real ratio_minus = 0.0;
    size_t particle_count = 0;
};

inline AphiGradDivSignEnergyMetrics hostGradDivSignEnergyMetrics(
    BaseParticles &particles, const AphiVariableNames &names, size_t total_real_particles, Real body_length,
    Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const auto in_core = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };

    AphiGradDivSignEnergyMetrics metrics;
    metrics.inner_a_plus_graddiv_a =
        hostAVectorVolDotProduct(particles, names.solution, names.v, total_real_particles, in_core);
    metrics.inner_a_minus_graddiv_a = -metrics.inner_a_plus_graddiv_a;
    metrics.div_a_norm_squared = hostDivANormSquared(particles, names, total_real_particles, in_core);
    if (metrics.div_a_norm_squared > TinyReal)
    {
        metrics.ratio_plus = metrics.inner_a_plus_graddiv_a / metrics.div_a_norm_squared;
        metrics.ratio_minus = metrics.inner_a_minus_graddiv_a / metrics.div_a_norm_squared;
    }
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_core(positions[i]))
        {
            metrics.particle_count += 1;
        }
    }
    return metrics;
}

struct AphiSourceCurrentDivergenceMetrics
{
    Real source_j_l2 = 0.0;
    Real source_div_j_l2 = 0.0;
    Real source_div_j_relative = 0.0;
    Real source_div_j_linf = 0.0;
    size_t source_region_count = 0;
    size_t total_particle_count = 0;
};

inline AphiSourceCurrentDivergenceMetrics hostSourceCurrentDivergenceMetrics(
    BaseParticles &particles, const AphiBlockNames &rhs_block, const std::string &div_j_real_name,
    const std::string &div_j_imag_name, const benchmark::AphiBoxRegion &source_region, size_t total_real_particles)
{
    syncVariableToHost<Vecd>(particles, rhs_block.a_real);
    syncVariableToHost<Vecd>(particles, rhs_block.a_imag);
    syncVariableToHost<Real>(particles, div_j_real_name);
    syncVariableToHost<Real>(particles, div_j_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");

    const Vecd *j_real = particles.getVariableDataByName<Vecd>(rhs_block.a_real);
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(rhs_block.a_imag);
    const Real *div_j_real = particles.getVariableDataByName<Real>(div_j_real_name);
    const Real *div_j_imag = particles.getVariableDataByName<Real>(div_j_imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    Real j_squared = 0.0;
    Real div_real_squared = 0.0;
    Real div_imag_squared = 0.0;
    Real grad_squared = 0.0;
    Real div_linf = 0.0;
    size_t source_count = 0;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        const bool in_source = benchmark::insideBoxRegion(positions[i], source_region);
        if (in_source)
        {
            source_count += 1;
        }
        if (!std::isfinite(div_j_real[i]) || !std::isfinite(div_j_imag[i]))
        {
            continue;
        }
        const Real vol_i = vol[i];
        j_squared += vol_i * (j_real[i].squaredNorm() + j_imag[i].squaredNorm());
        div_real_squared += vol_i * div_j_real[i] * div_j_real[i];
        div_imag_squared += vol_i * div_j_imag[i] * div_j_imag[i];
        grad_squared += vol_i * j_real[i].squaredNorm();
        div_linf = std::max(div_linf, std::abs(div_j_real[i]));
        div_linf = std::max(div_linf, std::abs(div_j_imag[i]));
    }

    AphiSourceCurrentDivergenceMetrics metrics;
    metrics.total_particle_count = total_real_particles;
    metrics.source_region_count = source_count;
    metrics.source_j_l2 = std::sqrt(j_squared);
    metrics.source_div_j_l2 = std::sqrt(div_real_squared + div_imag_squared);
    metrics.source_div_j_linf = div_linf;
    metrics.source_div_j_relative = metrics.source_div_j_l2 / (std::sqrt(grad_squared) + TinyReal);
    return metrics;
}

struct AphiCurlAManufacturedErrorMetrics
{
    Real b_error_l2 = 0.0;
    Real b_error_linf = 0.0;
    Real b_reference_l2 = 0.0;
    Real b_relative_error = 0.0;
    Real div_a_relative = 0.0;
    size_t particle_count = 0;
};

inline Real hostCoreDivARelative(BaseParticles &particles, const AphiVariableNames &names, size_t total_real_particles,
                                 Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_real);
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_imag);
    syncVariableToHost<Matd>(particles, names.solution.a_real + "Gradient");
    syncVariableToHost<Matd>(particles, names.solution.a_imag + "Gradient");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");

    const Real *div_a_real = particles.getVariableDataByName<Real>(names.diagnostic.div_a_real);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(names.diagnostic.div_a_imag);
    const Matd *grad_a_real = particles.getVariableDataByName<Matd>(names.solution.a_real + "Gradient");
    const Matd *grad_a_imag = particles.getVariableDataByName<Matd>(names.solution.a_imag + "Gradient");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    Real div_squared = 0.0;
    Real grad_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!std::isfinite(div_a_real[i]) || !grad_a_real[i].allFinite())
        {
            continue;
        }
        div_squared += vol[i] * (div_a_real[i] * div_a_real[i] + div_a_imag[i] * div_a_imag[i]);
        grad_squared += vol[i] * (grad_a_real[i].squaredNorm() + grad_a_imag[i].squaredNorm());
    }
    return std::sqrt(div_squared) / (std::sqrt(grad_squared) + TinyReal);
}

inline AphiCurlAManufacturedErrorMetrics hostCurlAManufacturedErrorMetrics(
    BaseParticles &particles, const AphiVariableNames &names, const std::string &b_real_name,
    const std::string &b_imag_name, const Vecd &reference_b_real, const Vecd &reference_b_imag,
    size_t total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Vecd>(particles, b_real_name);
    syncVariableToHost<Vecd>(particles, b_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");

    const Vecd *b_real = particles.getVariableDataByName<Vecd>(b_real_name);
    const Vecd *b_imag = particles.getVariableDataByName<Vecd>(b_imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    Real error_squared = 0.0;
    Real reference_squared = 0.0;
    Real error_linf = 0.0;
    size_t count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        count += 1;
        const Vecd err_re = b_real[i] - reference_b_real;
        const Vecd err_im = b_imag[i] - reference_b_imag;
        const Real vol_i = vol[i];
        error_squared += vol_i * (err_re.squaredNorm() + err_im.squaredNorm());
        reference_squared += vol_i * (reference_b_real.squaredNorm() + reference_b_imag.squaredNorm());
        error_linf = std::max(error_linf, err_re.norm());
        error_linf = std::max(error_linf, err_im.norm());
    }

    AphiCurlAManufacturedErrorMetrics metrics;
    metrics.particle_count = count;
    metrics.b_error_l2 = std::sqrt(error_squared);
    metrics.b_error_linf = error_linf;
    metrics.b_reference_l2 = std::sqrt(reference_squared);
    metrics.b_relative_error = metrics.b_error_l2 / (metrics.b_reference_l2 + TinyReal);
    metrics.div_a_relative =
        hostCoreDivARelative(particles, names, total_real_particles, body_length, body_height, body_width, core_shell);
    return metrics;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_A_GAUGE_DIAGNOSTIC_HELPERS_H
