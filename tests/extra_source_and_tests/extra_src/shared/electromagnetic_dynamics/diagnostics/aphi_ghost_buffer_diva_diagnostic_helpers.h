#ifndef APHI_GHOST_BUFFER_DIVA_DIAGNOSTIC_HELPERS_H
#define APHI_GHOST_BUFFER_DIVA_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_boundary_helpers.h"
#include "kernel_wendland_c2.h"

#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline bool isOutsidePhysicalDomain(const Vecd &position, Real body_length, Real body_height, Real body_width)
{
    return position[0] < 0.0 || position[0] > body_length || position[1] < 0.0 || position[1] > body_height ||
           position[2] < 0.0 || position[2] > body_width;
}

inline std::vector<Vecd> buildAnalyticGhostSamplePoints(Real dp, Real body_length, Real body_height, Real body_width,
                                                        Real ghost_layer_depth)
{
    std::vector<Vecd> ghost_points;
    const Real xmin = -ghost_layer_depth;
    const Real xmax = body_length + ghost_layer_depth;
    const Real ymin = -ghost_layer_depth;
    const Real ymax = body_height + ghost_layer_depth;
    const Real zmin = -ghost_layer_depth;
    const Real zmax = body_width + ghost_layer_depth;
    const size_t nx = static_cast<size_t>(std::ceil((xmax - xmin) / dp)) + 1;
    const size_t ny = static_cast<size_t>(std::ceil((ymax - ymin) / dp)) + 1;
    const size_t nz = static_cast<size_t>(std::ceil((zmax - zmin) / dp)) + 1;
    ghost_points.reserve(nx * ny * nz / 3);

    for (size_t ix = 0; ix <= nx; ++ix)
    {
        const Real x = xmin + ix * dp;
        for (size_t iy = 0; iy <= ny; ++iy)
        {
            const Real y = ymin + iy * dp;
            for (size_t iz = 0; iz <= nz; ++iz)
            {
                const Real z = zmin + iz * dp;
                const Vecd position(x, y, z);
                if (isOutsidePhysicalDomain(position, body_length, body_height, body_width))
                {
                    ghost_points.push_back(position);
                }
            }
        }
    }
    return ghost_points;
}

inline AphiDivAReductionMetrics hostRegionDivAGradDenMetricsFromStoredDivA(
    BaseParticles &particles, const std::string &div_a_real_name, const std::string &div_a_imag_name,
    const AphiVariableNames &names, const Vecd *positions, size_t total_real_particles, Real body_length,
    Real body_height, Real body_width, Real core_shell, const std::function<bool(const Vecd &)> &in_region)
{
    syncVariableToHost<Real>(particles, div_a_real_name);
    syncVariableToHost<Real>(particles, div_a_imag_name);
    syncVariableToHost<Matd>(particles, names.solution.a_real + "Gradient");
    syncVariableToHost<Matd>(particles, names.solution.a_imag + "Gradient");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Real *div_a_real = particles.getVariableDataByName<Real>(div_a_real_name);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(div_a_imag_name);
    const Matd *grad_a_real = particles.getVariableDataByName<Matd>(names.solution.a_real + "Gradient");
    const Matd *grad_a_imag = particles.getVariableDataByName<Matd>(names.solution.a_imag + "Gradient");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    Real div_real_squared = 0.0;
    Real div_imag_squared = 0.0;
    Real grad_squared = 0.0;
    Real div_linf = 0.0;
    size_t particle_count = 0;
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
        if (!isFiniteDivAGradEntry(div_a_real[i], grad_a_real[i]) ||
            !isFiniteDivAGradEntry(div_a_imag[i], grad_a_imag[i]))
        {
            continue;
        }
        particle_count += 1;
        const Real vol_i = vol[i];
        div_real_squared += vol_i * div_a_real[i] * div_a_real[i];
        div_imag_squared += vol_i * div_a_imag[i] * div_a_imag[i];
        grad_squared += vol_i * (grad_a_real[i].squaredNorm() + grad_a_imag[i].squaredNorm());
        div_linf = std::max(div_linf, std::abs(div_a_real[i]));
        div_linf = std::max(div_linf, std::abs(div_a_imag[i]));
    }
    return finalizeDivAReductionMetrics(div_real_squared, div_imag_squared, grad_squared, div_linf, particle_count);
}

struct AphiGhostBufferDivAMetrics
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real dp = 0.0;
    Real core_shell = 0.0;
    size_t ghost_point_count = 0;
    AphiDivAReductionMetrics baseline_boundary_grad_den{};
    AphiDivAReductionMetrics ghost_boundary_grad_den{};
    Real boundary_grad_den_reduction = 0.0;
};

inline Real hostSupplementalGhostDivAContribution(
    const Vecd &position_i, const Vecd &value_i, const Vecd *positions, size_t total_real_particles,
    const KernelWendlandC2 &kernel, const std::vector<Vecd> &ghost_positions,
    AphiDivFreeValidationFieldKind field_kind, Real ghost_vol, Real real_particle_merge_radius)
{
    Real divergence = 0.0;
    const auto exact_value = [&](const Vecd &position) {
        return divFreeValidationARealField(field_kind, position[0], position[1], position[2]);
    };
    for (const Vecd &ghost_position : ghost_positions)
    {
        bool covered_by_real = false;
        for (size_t j = 0; j != total_real_particles; ++j)
        {
            if ((ghost_position - positions[j]).norm() < real_particle_merge_radius)
            {
                covered_by_real = true;
                break;
            }
        }
        if (covered_by_real)
        {
            continue;
        }
        const Vecd displacement = position_i - ghost_position;
        if (displacement.squaredNorm() >= kernel.CutOffRadiusSqr())
        {
            continue;
        }
        const Real distance = displacement.norm();
        const Real dW_ij = kernel.dW(distance, displacement);
        const Vecd e_ij = displacement / (distance + TinyReal);
        const Vecd g_ij = -dW_ij * ghost_vol * e_ij;
        divergence += g_ij.dot(value_i - exact_value(ghost_position));
    }
    return divergence;
}

inline void hostApplySupplementalGhostDivACorrection(
    BaseParticles &particles, const AphiVariableNames &names, const Vecd *positions, size_t total_real_particles,
    Real body_length, Real body_height, Real body_width, Real dp, AphiDivFreeValidationFieldKind field_kind,
    const std::string &baseline_div_a_real, const std::string &baseline_div_a_imag,
    const std::string &ghost_div_a_real, const std::string &ghost_div_a_imag)
{
    syncAphiBlockToHost(particles, names.solution);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(names.solution.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(names.solution.a_imag);
    syncVariableToHost<Real>(particles, baseline_div_a_real);
    syncVariableToHost<Real>(particles, baseline_div_a_imag);

    const Real *baseline_re = particles.getVariableDataByName<Real>(baseline_div_a_real);
    const Real *baseline_im = particles.getVariableDataByName<Real>(baseline_div_a_imag);

    const Real reference_h = dp * 1.15;
    const KernelWendlandC2 kernel(reference_h);
    const std::vector<Vecd> ghost_positions =
        buildAnalyticGhostSamplePoints(dp, body_length, body_height, body_width, kernel.CutOffRadius());
    const Real ghost_vol = dp * dp * dp;
    const Real merge_radius = 0.5 * dp;

    particles.registerStateVariable<Real>(ghost_div_a_real, Real(0));
    particles.registerStateVariable<Real>(ghost_div_a_imag, Real(0));
    Real *ghost_re = particles.getVariableDataByName<Real>(ghost_div_a_real);
    Real *ghost_im = particles.getVariableDataByName<Real>(ghost_div_a_imag);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        ghost_re[i] =
            baseline_re[i] + hostSupplementalGhostDivAContribution(positions[i], a_real[i], positions, total_real_particles,
                                                                   kernel, ghost_positions, field_kind, ghost_vol, merge_radius);
        ghost_im[i] =
            baseline_im[i] + hostSupplementalGhostDivAContribution(positions[i], a_imag[i], positions, total_real_particles,
                                                                   kernel, ghost_positions, field_kind, ghost_vol, merge_radius);
    }
}

inline AphiGhostBufferDivAMetrics runGhostBufferDivAMetrics(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, AphiDivFreeValidationFieldKind field_kind)
{
    const Real boundary_width = 3.0 * dp_0;
    AphiGhostBufferDivAMetrics metrics;
    metrics.field_kind = field_kind;
    metrics.dp = dp_0;
    metrics.core_shell = core_shell;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivFreeValidationAphiFieldCK> assign_field(
        test_body.body, names.solution, field_kind);

    initialize_aphi.exec();
    set_material.exec();
    assign_field.exec();
    test_body.updateRelations();

    execBodyVectorAGradientForDivADiagnostic(test_body.body, test_body.inner(), names);
    execBodyPairwiseDivADiagnosticPipeline(test_body.body, test_body.inner(), names);

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();
    const auto in_boundary = [&](const Vecd &position) {
        return !isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };

    const KernelWendlandC2 kernel(dp_0 * 1.15);
    metrics.ghost_point_count = buildAnalyticGhostSamplePoints(dp_0, body_length, body_height, body_width,
                                                             kernel.CutOffRadius())
                                    .size();

    metrics.baseline_boundary_grad_den = hostRegionDivAGradDenMetricsFromStoredDivA(
        particles, names.diagnostic.div_a_real, names.diagnostic.div_a_imag, names, positions, total_real_particles,
        body_length, body_height, body_width, core_shell, in_boundary);

    hostApplySupplementalGhostDivACorrection(particles, names, positions, total_real_particles, body_length, body_height,
                                           body_width, dp_0, field_kind, names.diagnostic.div_a_real,
                                           names.diagnostic.div_a_imag, "DivARealGhostHost", "DivAImagGhostHost");
    metrics.ghost_boundary_grad_den = hostRegionDivAGradDenMetricsFromStoredDivA(
        particles, "DivARealGhostHost", "DivAImagGhostHost", names, positions, total_real_particles, body_length,
        body_height, body_width, core_shell, in_boundary);

    const Real baseline = metrics.baseline_boundary_grad_den.div_a_relative;
    const Real ghost = metrics.ghost_boundary_grad_den.div_a_relative;
    metrics.boundary_grad_den_reduction = (baseline - ghost) / (baseline + TinyReal);
    return metrics;
}

inline void printGhostBufferDivAMetrics(const char *test_name, const AphiGhostBufferDivAMetrics &metrics)
{
    std::cout << test_name << " field=" << divFreeValidationFieldName(metrics.field_kind) << " dp=" << metrics.dp
              << " core_shell=" << metrics.core_shell << " ghost_points=" << metrics.ghost_point_count
              << " baseline_boundary_gradDen_rel=" << metrics.baseline_boundary_grad_den.div_a_relative
              << " ghost_boundary_gradDen_rel=" << metrics.ghost_boundary_grad_den.div_a_relative
              << " boundary_reduction=" << metrics.boundary_grad_den_reduction << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GHOST_BUFFER_DIVA_DIAGNOSTIC_HELPERS_H
