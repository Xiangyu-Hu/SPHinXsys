#ifndef APHI_PAIRWISE_DIVA_BOUNDARY_HELPERS_H
#define APHI_PAIRWISE_DIVA_BOUNDARY_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_pairwise_diva_root_cause_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiPairwiseDivABoundaryMetrics
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real dp = 0.0;
    Real boundary_width = 0.0;
    Real core_shell = 0.0;
    AphiDivAReductionMetrics core_pairwise_grad_den{};
    AphiDivAReductionMetrics boundary_pairwise_grad_den{};
    AphiDivAReductionMetrics core_b_corrected_grad_den{};
    AphiDivAReductionMetrics boundary_b_corrected_grad_den{};
    Real boundary_div_energy_fraction = 0.0;
    Real mean_boundary_distance_to_wall = 0.0;
    Real mean_boundary_div_a_abs = 0.0;
    size_t core_particles = 0;
    size_t boundary_particles = 0;
};

inline AphiPairwiseDivABoundaryMetrics runPairwiseDivABoundaryMetrics(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, AphiDivFreeValidationFieldKind field_kind)
{
    AphiPairwiseDivABoundaryMetrics metrics;
    metrics.field_kind = field_kind;
    metrics.dp = dp_0;
    metrics.boundary_width = boundary_width;
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
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_real);
    syncVariableToHost<Real>(particles, names.diagnostic.div_a_imag);
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *div_a_real = particles.getVariableDataByName<Real>(names.diagnostic.div_a_real);
    const Real *div_a_imag = particles.getVariableDataByName<Real>(names.diagnostic.div_a_imag);
    const size_t total_real_particles = particles.TotalRealParticles();

    const auto in_core = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    const auto in_boundary = [&](const Vecd &position) {
        return !isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };

    metrics.core_pairwise_grad_den = hostBodyPairwiseDivWithGradDenominatorMetrics(
        particles, names, positions, total_real_particles, body_length, body_height, body_width, core_shell, in_core);
    metrics.boundary_pairwise_grad_den = hostBodyPairwiseDivWithGradDenominatorMetrics(
        particles, names, positions, total_real_particles, body_length, body_height, body_width, core_shell, in_boundary);
    metrics.core_particles = metrics.core_pairwise_grad_den.particle_count;
    metrics.boundary_particles = metrics.boundary_pairwise_grad_den.particle_count;

    const AphiDivAReductionMetrics global_pairwise = hostBodyDivAReductionMetrics(
        particles, names, positions, total_real_particles, {}, AphiDivADiagnosticMode::PairwiseUncorrected);
    const Real total_div_energy = global_pairwise.div_a_L2 * global_pairwise.div_a_L2;
    const Real boundary_div_energy = metrics.boundary_pairwise_grad_den.div_a_L2 * metrics.boundary_pairwise_grad_den.div_a_L2;
    metrics.boundary_div_energy_fraction = boundary_div_energy / (total_div_energy + TinyReal);

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names, AphiDivADiagnosticMode::BCorrectedTrace);
    metrics.core_b_corrected_grad_den =
        hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, in_core,
                                     AphiDivADiagnosticMode::BCorrectedTrace);
    metrics.boundary_b_corrected_grad_den =
        hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, in_boundary,
                                     AphiDivADiagnosticMode::BCorrectedTrace);

    Real distance_sum = 0.0;
    Real div_abs_sum = 0.0;
    size_t boundary_count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (in_core(positions[i]))
        {
            continue;
        }
        boundary_count += 1;
        distance_sum += distanceToBoundary(positions[i], body_length, body_height, body_width);
        div_abs_sum += std::abs(div_a_real[i]) + std::abs(div_a_imag[i]);
    }
    metrics.mean_boundary_distance_to_wall = boundary_count > 0 ? distance_sum / boundary_count : 0.0;
    metrics.mean_boundary_div_a_abs = boundary_count > 0 ? div_abs_sum / boundary_count : 0.0;
    return metrics;
}

inline void printPairwiseDivABoundaryMetrics(const char *test_name, const AphiPairwiseDivABoundaryMetrics &metrics)
{
    std::cout << test_name << " field=" << divFreeValidationFieldName(metrics.field_kind) << " dp=" << metrics.dp
              << " boundary_width=" << metrics.boundary_width << " core_shell=" << metrics.core_shell
              << " core_N=" << metrics.core_particles << " boundary_N=" << metrics.boundary_particles
              << " core_pairwise_gradDen_rel=" << metrics.core_pairwise_grad_den.div_a_relative
              << " boundary_pairwise_gradDen_rel=" << metrics.boundary_pairwise_grad_den.div_a_relative
              << " core_B_gradDen_rel=" << metrics.core_b_corrected_grad_den.div_a_relative
              << " boundary_B_gradDen_rel=" << metrics.boundary_b_corrected_grad_den.div_a_relative
              << " boundary_div_energy_frac=" << metrics.boundary_div_energy_fraction
              << " mean_boundary_dist=" << metrics.mean_boundary_distance_to_wall
              << " mean_boundary_divA_abs=" << metrics.mean_boundary_div_a_abs << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PAIRWISE_DIVA_BOUNDARY_HELPERS_H
