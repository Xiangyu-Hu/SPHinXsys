#ifndef APHI_INNER_A_DIVERGENCE_PENALTY_ETA_SWEEP_HELPERS_H
#define APHI_INNER_A_DIVERGENCE_PENALTY_ETA_SWEEP_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_sweep_helpers.h"

#include <algorithm>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiCoreOperatorScaleMetrics
{
    Real median_abs_laplace_a_diag = 0.0;
    Real median_graddiv_block_norm = 0.0;
    size_t core_particles = 0;
};

struct AphiInnerADivergencePenaltyEtaSweepRow : public AphiInnerADivergencePenaltySweepRow
{
    Real eta_a = 0.0;
    Real penalty_to_laplace_diag_ratio = 0.0;
};

inline Real hostCoreMedianAbsReal(BaseParticles &particles, const std::string &field_name, const Vecd *positions,
                                  size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                                  Real core_shell)
{
    syncVariableToHost<Real>(particles, field_name);
    const Real *values = particles.getVariableDataByName<Real>(field_name);
    std::vector<Real> core_values;
    core_values.reserve(total_real_particles);
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
        core_values.push_back(std::abs(values[i]));
    }
    if (core_values.empty())
    {
        return 0.0;
    }
    std::sort(core_values.begin(), core_values.end());
    return core_values[core_values.size() / 2];
}

inline Real hostCoreMedianMatNorm(BaseParticles &particles, const std::string &field_name, const Vecd *positions,
                                  size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                                  Real core_shell)
{
    syncVariableToHost<Matd>(particles, field_name);
    const Matd *values = particles.getVariableDataByName<Matd>(field_name);
    std::vector<Real> core_values;
    core_values.reserve(total_real_particles);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!values[i].allFinite())
        {
            continue;
        }
        core_values.push_back(values[i].norm());
    }
    if (core_values.empty())
    {
        return 0.0;
    }
    std::sort(core_values.begin(), core_values.end());
    return core_values[core_values.size() / 2];
}

inline AphiCoreOperatorScaleMetrics hostCoreOperatorScaleMetrics(
    AphiLhsTestBody &test_body, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, Real omega, Real phi_gauge_penalty, int ac, char *av[])
{
    (void)ac;
    (void)av;
    AphiVariableNames names;
    const AphiBlockJacobiDiagonalNames diag_names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = true;
    options.a_divergence_penalty = 1.0;

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi_diagonal(
        test_body.inner(), names.material, omega, options, diag_names);

    initialize_aphi.exec();
    set_material.exec();
    test_body.updateRelations();
    compute_jacobi_diagonal.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    AphiCoreOperatorScaleMetrics metrics;
    metrics.median_abs_laplace_a_diag =
        hostCoreMedianAbsReal(particles, diag_names.laplace_a_diag, positions, total_real_particles, body_length,
                              body_height, body_width, core_shell);
    metrics.median_graddiv_block_norm =
        hostCoreMedianMatNorm(particles, diag_names.graddiv_a_block, positions, total_real_particles, body_length,
                              body_height, body_width, core_shell);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            metrics.core_particles += 1;
        }
    }
    return metrics;
}

inline Real lambdaAFromEtaA(Real eta_a, const AphiCoreOperatorScaleMetrics &scale_metrics)
{
    // eta_a is penalty_to_laplace_diag_ratio; see AphiADivergencePenaltyResearchDefaults.
    return eta_a * scale_metrics.median_abs_laplace_a_diag / (scale_metrics.median_graddiv_block_norm + TinyReal);
}

inline AphiInnerADivergencePenaltyEtaSweepRow toEtaSweepRow(const AphiInnerADivergencePenaltySweepRow &base)
{
    AphiInnerADivergencePenaltyEtaSweepRow row;
    static_cast<AphiInnerADivergencePenaltySweepRow &>(row) = base;
    return row;
}

inline AphiInnerADivergencePenaltyEtaSweepRow runInnerImpressedEtaSweepRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, const benchmark::AphiBoxRegion &source_region,
    const Vecd &current_real, const Vecd &current_imag, Real impressed_current_amplitude,
    UnsignedInt restart_dimension = 50, UnsignedInt max_outer_iterations = 100)
{
    AphiInnerADivergencePenaltyEtaSweepRow row;
    row.eta_a = eta_a;
    if (eta_a <= TinyReal)
    {
        row = toEtaSweepRow(runInnerADivergencePenaltySweepRow(ac, av, dp_0, body_length, body_height, body_width,
                                                              boundary_width, core_shell, sigma, nu, omega,
                                                              phi_gauge_penalty, 0.0, false, tolerance, source_region,
                                                              current_real, current_imag, impressed_current_amplitude,
                                                              restart_dimension, max_outer_iterations));
        row.eta_a = eta_a;
        row.penalty_to_laplace_diag_ratio = 0.0;
        return row;
    }

    const Real lambda_a = lambdaAFromEtaA(eta_a, scale_metrics);
    row = toEtaSweepRow(runInnerADivergencePenaltySweepRow(ac, av, dp_0, body_length, body_height, body_width,
                                                           boundary_width, core_shell, sigma, nu, omega,
                                                           phi_gauge_penalty, lambda_a, true, tolerance, source_region,
                                                           current_real, current_imag, impressed_current_amplitude,
                                                           restart_dimension, max_outer_iterations));
    row.eta_a = eta_a;
    row.a_divergence_penalty = lambda_a;
    row.penalty_to_laplace_diag_ratio = eta_a;
    return row;
}

inline AphiInnerADivergencePenaltyEtaSweepRow runInnerSolenoidalEtaSweepRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, const benchmark::AphiBoxRegion &source_region,
    Real solenoidal_current_amplitude, Real imag_to_real_ratio = 0.1, UnsignedInt restart_dimension = 50,
    UnsignedInt max_outer_iterations = 100)
{
    AphiInnerADivergencePenaltyEtaSweepRow row;
    row.eta_a = eta_a;
    if (eta_a <= TinyReal)
    {
        row = toEtaSweepRow(runInnerSolenoidalADivergencePenaltySweepRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
            phi_gauge_penalty, 0.0, false, tolerance, source_region, solenoidal_current_amplitude, imag_to_real_ratio,
            restart_dimension, max_outer_iterations));
        row.eta_a = eta_a;
        row.penalty_to_laplace_diag_ratio = 0.0;
        return row;
    }

    const Real lambda_a = lambdaAFromEtaA(eta_a, scale_metrics);
    row = toEtaSweepRow(runInnerSolenoidalADivergencePenaltySweepRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, lambda_a, true, tolerance, source_region, solenoidal_current_amplitude, imag_to_real_ratio,
        restart_dimension, max_outer_iterations));
    row.eta_a = eta_a;
    row.a_divergence_penalty = lambda_a;
    row.penalty_to_laplace_diag_ratio = eta_a;
    return row;
}

inline void printInnerADivergencePenaltyEtaSweepRow(const char *test_name,
                                                    const AphiInnerADivergencePenaltyEtaSweepRow &row)
{
    std::cout << test_name << " eta_A=" << row.eta_a << " lambda_A=" << row.a_divergence_penalty
              << " penalty_to_laplace_diag_ratio=" << row.penalty_to_laplace_diag_ratio
              << " use_a_penalty=" << (row.use_a_divergence_penalty ? 1 : 0) << " converged=" << (row.converged ? 1 : 0)
              << " outer=" << row.outer_iterations << " true_rel=" << row.final_true_relative_residual
              << " res_A_frac=" << row.final_res_a_fraction << " res_phi_frac=" << row.final_res_phi_fraction
              << " div_A_relative=" << row.global_div_a.div_a_relative
              << " source_div_A_relative=" << row.source_region_div_a.div_a_relative
              << " div_A_level=" << divAGaugeDiagnosticLevel(row.global_div_a.div_a_relative) << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_INNER_A_DIVERGENCE_PENALTY_ETA_SWEEP_HELPERS_H
