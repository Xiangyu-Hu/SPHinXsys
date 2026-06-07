#ifndef APHI_PAIRWISE_DIVA_ROOT_CAUSE_HELPERS_H
#define APHI_PAIRWISE_DIVA_ROOT_CAUSE_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline void execBodyVectorAGradientForDivADiagnostic(SPHBody &body, Inner<> &inner, const AphiVariableNames &names)
{
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(inner);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient(inner, names.solution.a_real);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient(inner, names.solution.a_imag);
    linear_correction_matrix.exec();
    a_real_gradient.exec();
    a_imag_gradient.exec();
}

inline AphiDivAReductionMetrics hostBodyPairwiseDivWithGradDenominatorMetrics(
    BaseParticles &particles, const AphiVariableNames &names, const Vecd *positions, size_t total_real_particles,
    Real body_length, Real body_height, Real body_width, Real core_shell,
    const std::function<bool(const Vecd &)> &in_region = {})
{
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_real);
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_imag);
    syncVariableToHost<Matd>(particles, names.solution.a_real + "Gradient");
    syncVariableToHost<Matd>(particles, names.solution.a_imag + "Gradient");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Real *div_a_real = particles.getVariableDataByName<Real>(names.diagnostic.div_a_real);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(names.diagnostic.div_a_imag);
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
        if (in_region)
        {
            if (!in_region(positions[i]))
            {
                continue;
            }
        }
        else if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
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

struct AphiPairwiseDivARootCauseMetrics
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real dp = 0.0;
    Real core_shell = 0.0;
    AphiDivAReductionMetrics pairwise_a_norm{};
    AphiDivAReductionMetrics pairwise_grad_den{};
    AphiDivAReductionMetrics b_corrected_grad_norm{};
    AphiDivAReductionMetrics core_pairwise_a_norm{};
    AphiDivAReductionMetrics boundary_pairwise_a_norm{};
    Real boundary_div_energy_fraction = 0.0;
    size_t core_particles = 0;
    size_t boundary_particles = 0;
};

inline AphiPairwiseDivARootCauseMetrics runPairwiseDivARootCauseMetrics(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, AphiDivFreeValidationFieldKind field_kind)
{
    const Real boundary_width = 3.0 * dp_0;
    AphiPairwiseDivARootCauseMetrics metrics;
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

    metrics.pairwise_a_norm = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, {},
                                                           AphiDivADiagnosticMode::PairwiseUncorrected);
    metrics.pairwise_grad_den = hostBodyPairwiseDivWithGradDenominatorMetrics(
        particles, names, positions, total_real_particles, body_length, body_height, body_width, core_shell, {});

    const auto in_core = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    const auto in_boundary = [&](const Vecd &position) {
        return !isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    metrics.core_pairwise_a_norm =
        hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, in_core,
                                     AphiDivADiagnosticMode::PairwiseUncorrected);
    metrics.boundary_pairwise_a_norm =
        hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, in_boundary,
                                     AphiDivADiagnosticMode::PairwiseUncorrected);
    metrics.core_particles = metrics.core_pairwise_a_norm.particle_count;
    metrics.boundary_particles = metrics.boundary_pairwise_a_norm.particle_count;

    const Real total_div_energy = metrics.pairwise_a_norm.div_a_L2 * metrics.pairwise_a_norm.div_a_L2;
    const Real boundary_div_energy = metrics.boundary_pairwise_a_norm.div_a_L2 * metrics.boundary_pairwise_a_norm.div_a_L2;
    metrics.boundary_div_energy_fraction = boundary_div_energy / (total_div_energy + TinyReal);

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names, AphiDivADiagnosticMode::BCorrectedTrace);
    metrics.b_corrected_grad_norm = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, {},
                                                               AphiDivADiagnosticMode::BCorrectedTrace);
    return metrics;
}

inline void printPairwiseDivARootCauseMetrics(const char *test_name, const AphiPairwiseDivARootCauseMetrics &metrics)
{
    std::cout << test_name << " field=" << divFreeValidationFieldName(metrics.field_kind) << " dp=" << metrics.dp
              << " core_shell=" << metrics.core_shell << " core_N=" << metrics.core_particles
              << " boundary_N=" << metrics.boundary_particles << " pairwise_divA_rel_Anorm="
              << metrics.pairwise_a_norm.div_a_relative
              << " pairwise_divA_rel_gradDen=" << metrics.pairwise_grad_den.div_a_relative
              << " B_divA_rel_gradDen=" << metrics.b_corrected_grad_norm.div_a_relative
              << " core_pairwise_rel=" << metrics.core_pairwise_a_norm.div_a_relative
              << " boundary_pairwise_rel=" << metrics.boundary_pairwise_a_norm.div_a_relative
              << " boundary_div_energy_frac=" << metrics.boundary_div_energy_fraction << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PAIRWISE_DIVA_ROOT_CAUSE_HELPERS_H
