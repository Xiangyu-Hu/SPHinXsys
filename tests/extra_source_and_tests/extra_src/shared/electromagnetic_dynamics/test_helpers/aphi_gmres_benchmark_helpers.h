#ifndef APHI_GMRES_BENCHMARK_HELPERS_H
#define APHI_GMRES_BENCHMARK_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_gmres_robustness_sweep_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"
#include "electromagnetic_dynamics/benchmark/aphi_team7_canonical_case_ck.h"
#include "electromagnetic_dynamics/benchmark/aphi_cold_crucible_case_ck.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <ostream>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiGMRESBenchmarkRunResult
{
    AphiManufacturedRobustnessResult solver;
    Real core_relative_solution_error = 0.0;
    Real source_region_solution_block_max = 0.0;
    Real exterior_solution_block_max = 0.0;
    Real conductor_solution_block_max = 0.0;
    Real coil_solution_block_max = 0.0;
};

struct AphiComponentWiseNorms
{
    Real a_real_l2 = 0.0;
    Real a_imag_l2 = 0.0;
    Real phi_real_l2 = 0.0;
    Real phi_imag_l2 = 0.0;
};

enum class AphiBenchmarkMaterialRegion
{
    AllCore,
    Conductor,
    Coil,
    Air
};

inline bool team7ParticleInRegion(const Vecd &position, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                  AphiBenchmarkMaterialRegion region)
{
    switch (region)
    {
    case AphiBenchmarkMaterialRegion::Conductor:
        return benchmark::insideBoxRegion(position, layout.conductor);
    case AphiBenchmarkMaterialRegion::Coil:
        return benchmark::insideBoxRegion(position, layout.coil);
    case AphiBenchmarkMaterialRegion::Air:
        return !benchmark::insideBoxRegion(position, layout.conductor) &&
               !benchmark::insideBoxRegion(position, layout.coil);
    case AphiBenchmarkMaterialRegion::AllCore:
    default:
        return true;
    }
}

inline AphiComponentWiseNorms hostTrueResidualComponentNorms(
    BaseParticles &particles, const AphiBlockNames &residual_block, const Vecd *positions, size_t total_real_particles,
    Real body_length, Real body_height, Real body_width, Real core_shell,
    const std::function<bool(const Vecd &)> &include_particle)
{
    syncAphiBlockToHost(particles, residual_block);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *res_a_real = particles.getVariableDataByName<Vecd>(residual_block.a_real);
    const Vecd *res_a_imag = particles.getVariableDataByName<Vecd>(residual_block.a_imag);
    const Real *res_phi_real = particles.getVariableDataByName<Real>(residual_block.phi_real);
    const Real *res_phi_imag = particles.getVariableDataByName<Real>(residual_block.phi_imag);

    Real sum_a_real = 0.0;
    Real sum_a_imag = 0.0;
    Real sum_phi_real = 0.0;
    Real sum_phi_imag = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!include_particle(positions[i]))
        {
            continue;
        }
        sum_a_real += vol[i] * res_a_real[i].squaredNorm();
        sum_a_imag += vol[i] * res_a_imag[i].squaredNorm();
        sum_phi_real += vol[i] * res_phi_real[i] * res_phi_real[i];
        sum_phi_imag += vol[i] * res_phi_imag[i] * res_phi_imag[i];
    }

    AphiComponentWiseNorms norms;
    norms.a_real_l2 = std::sqrt(sum_a_real);
    norms.a_imag_l2 = std::sqrt(sum_a_imag);
    norms.phi_real_l2 = std::sqrt(sum_phi_real);
    norms.phi_imag_l2 = std::sqrt(sum_phi_imag);
    return norms;
}

inline Real hostTrueResidualBlockNormInRegion(
    BaseParticles &particles, const AphiBlockNames &residual_block, const Vecd *positions, size_t total_real_particles,
    Real body_length, Real body_height, Real body_width, Real core_shell,
    const std::function<bool(const Vecd &)> &include_particle)
{
    const AphiComponentWiseNorms norms =
        hostTrueResidualComponentNorms(particles, residual_block, positions, total_real_particles, body_length,
                                       body_height, body_width, core_shell, include_particle);
    return std::sqrt(norms.a_real_l2 * norms.a_real_l2 + norms.a_imag_l2 * norms.a_imag_l2 +
                     norms.phi_real_l2 * norms.phi_real_l2 + norms.phi_imag_l2 * norms.phi_imag_l2);
}

inline void logTrueResidualDiagnostics(std::ostream &out, BaseParticles &particles, const AphiVariableNames &names,
                                       const Vecd *positions, size_t total_real_particles, Real body_length,
                                       Real body_height, Real body_width, Real core_shell, Real reference_norm,
                                       const benchmark::AphiTeam7LikeUnitBoxLayout &layout, const char *label)
{
    const auto all_core = [](const Vecd &) { return true; };
    const AphiComponentWiseNorms component_norms = hostTrueResidualComponentNorms(
        particles, names.true_residual, positions, total_real_particles, body_length, body_height, body_width,
        core_shell, all_core);

    out << " " << label << "_component_rel"
        << " a_real=" << (component_norms.a_real_l2 / (reference_norm + TinyReal))
        << " a_imag=" << (component_norms.a_imag_l2 / (reference_norm + TinyReal))
        << " phi_real=" << (component_norms.phi_real_l2 / (reference_norm + TinyReal))
        << " phi_imag=" << (component_norms.phi_imag_l2 / (reference_norm + TinyReal));

    const auto region_norm = [&](AphiBenchmarkMaterialRegion region) {
        return hostTrueResidualBlockNormInRegion(
            particles, names.true_residual, positions, total_real_particles, body_length, body_height, body_width,
            core_shell, [&](const Vecd &position) { return team7ParticleInRegion(position, layout, region); });
    };

    out << " " << label << "_region_rel"
        << " conductor=" << (region_norm(AphiBenchmarkMaterialRegion::Conductor) / (reference_norm + TinyReal))
        << " coil=" << (region_norm(AphiBenchmarkMaterialRegion::Coil) / (reference_norm + TinyReal))
        << " air=" << (region_norm(AphiBenchmarkMaterialRegion::Air) / (reference_norm + TinyReal));
}

/** TEAM7-like GMRES run with shared max workspace (supports m <= workspace dimension). */
inline AphiGMRESResult runTeam7LikeGMRESDirect(
    AphiLhsTestBody &test_body, const AphiVariableNames &names, const AphiGMRESWorkspaceNames &workspace,
    const AphiLhsAssemblyOptions &options, StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> &zero_solution,
    Real tolerance, UnsignedInt restart_dimension, UnsignedInt max_outer_iterations, bool use_block_jacobi_preconditioner,
    AphiBlockJacobiPreconditionerKind block_jacobi_kind = AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8)
{
    zero_solution.exec();
    test_body.updateRelations();

    AphiGMRESSolverOptions gmres_options =
        defaultGMRESConvergenceOptions(tolerance, restart_dimension, max_outer_iterations);
    gmres_options.use_block_jacobi_preconditioner = use_block_jacobi_preconditioner;
    gmres_options.block_jacobi_kind = block_jacobi_kind;
    AphiGMRESSolverCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, workspace, options,
                                                  gmres_options);
    return solver.solve();
}

inline Real coreBlockMaxAbsNormInRegion(BaseParticles &particles, const AphiBlockNames &block_names,
                                        const Vecd *positions, size_t total_real_particles, Real body_length,
                                        Real body_height, Real body_width, Real core_shell,
                                        const benchmark::AphiBoxRegion &region)
{
    syncAphiBlockToHost(particles, block_names);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block_names.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block_names.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block_names.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block_names.phi_imag);
    Real max_norm = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!benchmark::insideBoxRegion(positions[i], region))
        {
            continue;
        }
        max_norm = std::max(max_norm, a_real[i].norm());
        max_norm = std::max(max_norm, a_imag[i].norm());
        max_norm = std::max(max_norm, std::abs(phi_real[i]));
        max_norm = std::max(max_norm, std::abs(phi_imag[i]));
    }
    return max_norm;
}

inline Real coreBlockMaxAbsNormInRegionPredicate(BaseParticles &particles, const AphiBlockNames &block_names,
                                                 const Vecd *positions, size_t total_real_particles, Real body_length,
                                                 Real body_height, Real body_width, Real core_shell,
                                                 const std::function<bool(const Vecd &)> &in_region)
{
    syncAphiBlockToHost(particles, block_names);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block_names.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block_names.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block_names.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block_names.phi_imag);
    Real max_norm = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!in_region(positions[i]))
        {
            continue;
        }
        max_norm = std::max(max_norm, a_real[i].norm());
        max_norm = std::max(max_norm, a_imag[i].norm());
        max_norm = std::max(max_norm, std::abs(phi_real[i]));
        max_norm = std::max(max_norm, std::abs(phi_imag[i]));
    }
    return max_norm;
}

inline Real coreBlockMaxAbsNormOutsideRegion(BaseParticles &particles, const AphiBlockNames &block_names,
                                             const Vecd *positions, size_t total_real_particles, Real body_length,
                                             Real body_height, Real body_width, Real core_shell,
                                             const benchmark::AphiBoxRegion &region)
{
    syncAphiBlockToHost(particles, block_names);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block_names.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block_names.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block_names.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block_names.phi_imag);
    Real max_norm = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (benchmark::insideBoxRegion(positions[i], region))
        {
            continue;
        }
        max_norm = std::max(max_norm, a_real[i].norm());
        max_norm = std::max(max_norm, a_imag[i].norm());
        max_norm = std::max(max_norm, std::abs(phi_real[i]));
        max_norm = std::max(max_norm, std::abs(phi_imag[i]));
    }
    return max_norm;
}

inline Real coreRelativeBlockError(BaseParticles &particles, const AphiBlockNames &approx_block,
                                   const AphiBlockNames &exact_block, const Vecd *positions,
                                   size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                                   Real core_shell)
{
    const Real error_norm = coreBlockMaxAbsDifference(particles, approx_block, exact_block, positions,
                                                      total_real_particles, body_length, body_height, body_width,
                                                      core_shell);
    const Real exact_norm = coreBlockOperatorNorm(particles, exact_block, positions, total_real_particles, body_length,
                                                  body_height, body_width, core_shell);
    return error_norm / (exact_norm + TinyReal);
}

/** ||K(u_approx - u_exact)|| / ||rhs|| on core particles — discrete MMS defect (gauge-aware). */
inline Real coreDiscreteMmsOperatorDefect(AphiLhsTestBody &test_body, const AphiVariableNames &names,
                                          const AphiLhsAssemblyOptions &options,
                                          const AphiBlockNames &approx_block, const AphiBlockNames &exact_block)
{
    StateDynamics<MainExecutionPolicy, AphiBlockLinearCombinationCK> form_difference(
        test_body.body, names.t, Real(1), Real(-1), approx_block, exact_block);
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_difference(test_body.body, test_body.inner(), names.t, names.lhs,
                                                                  names.material, options.omega, options);

    form_difference.exec();
    apply_difference.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Real defect_norm = hostBlockNorm(particles, names.lhs, total_real_particles);
    const Real rhs_norm_l2 = hostBlockNorm(particles, names.rhs, total_real_particles);
    return defect_norm / (rhs_norm_l2 + TinyReal);
}

/** Vol-weighted true-residual L2 norm on core particles within |x - x_interface| <= band_half_width. */
inline Real hostInterfaceBandTrueResidualNorm(BaseParticles &particles, const AphiBlockNames &residual_block,
                                              const Vecd *positions, size_t total_real_particles, Real body_length,
                                              Real body_height, Real body_width, Real core_shell, Real x_interface,
                                              Real band_half_width)
{
    const auto in_band = [&](const Vecd &position) {
        return std::abs(position[0] - x_interface) <= band_half_width;
    };
    return hostTrueResidualBlockNormInRegion(particles, residual_block, positions, total_real_particles, body_length,
                                             body_height, body_width, core_shell, in_band);
}

/** Frequency continuation: ramp omega, warm-starting from the previous solution. */
inline AphiMatrixFreeSolverResult runGMRESWithOmegaContinuation(
    AphiLhsTestBody &test_body, const AphiVariableNames &names, AphiLhsAssemblyOptions &options,
    const std::vector<Real> &omega_steps, Real tolerance, UnsignedInt restart_dimension = 30,
    UnsignedInt max_outer_iterations = 25, bool use_block_jacobi_preconditioner = true,
    bool require_step_convergence = true)
{
    AphiMatrixFreeSolverResult last_result;
    for (const Real omega_step : omega_steps)
    {
        options.omega = omega_step;
        AphiMatrixFreeSolverOptions solver_options =
            defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
        solver_options.gmres.use_block_jacobi_preconditioner = use_block_jacobi_preconditioner;
        AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                         solver_options);
        test_body.updateRelations();
        last_result = solver.solve();
        if (require_step_convergence && !gmresConvergencePassed(last_result, tolerance))
        {
            return last_result;
        }
    }
    return last_result;
}

/** Ramp conductor sigma from start to target; warm-start between steps. */
inline AphiMatrixFreeSolverResult runGMRESWithConductorSigmaContinuation(
    AphiLhsTestBody &test_body, const AphiVariableNames &names, AphiLhsAssemblyOptions &options,
    benchmark::AphiTeam7LikeUnitBoxLayout &layout, Real tolerance, const std::vector<Real> &conductor_sigma_steps,
    UnsignedInt restart_dimension = 30, UnsignedInt max_outer_iterations = 40,
    bool use_block_jacobi_preconditioner = true)
{
    AphiMatrixFreeSolverResult last_result;
    for (const Real sigma_step : conductor_sigma_steps)
    {
        layout.conductor_material.sigma = sigma_step;
        StateDynamics<MainExecutionPolicy, benchmark::AssignTeam7LikeRegionMaterialsCK> assign_material(
            test_body.body, layout, names.material);
        assign_material.exec();
        test_body.updateRelations();

        AphiMatrixFreeSolverOptions solver_options =
            defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
        solver_options.gmres.use_block_jacobi_preconditioner = use_block_jacobi_preconditioner;
        AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                         solver_options);
        last_result = solver.solve();
        if (!gmresConvergencePassed(last_result, tolerance))
        {
            return last_result;
        }
    }
    return last_result;
}

inline Real hostVolWeightedSum(BaseParticles &particles, const std::string &variable_name, size_t total_real_particles,
                               const std::function<bool(const Vecd &)> &include_particle,
                               const Vecd *positions = nullptr)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    if (positions != nullptr)
    {
        syncVariableToHost<Vecd>(particles, "Position");
    }
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *position = positions != nullptr ? particles.getVariableDataByName<Vecd>("Position") : nullptr;
    Real weighted_sum = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (position != nullptr && !include_particle(position[i]))
        {
            continue;
        }
        weighted_sum += vol[i] * values[i];
    }
    return weighted_sum;
}

inline Real hostTeam7RegionJoulePower(BaseParticles &particles, const Vecd *positions, size_t total_real_particles,
                                      const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                      AphiBenchmarkMaterialRegion region, const std::string &joule_source_name)
{
    return hostVolWeightedSum(
        particles, joule_source_name, total_real_particles,
        [&](const Vecd &position) { return team7ParticleInRegion(position, layout, region); }, positions);
}

struct AphiJoulePowerMetrics
{
    Real total = 0.0;
    Real conductor = 0.0;
    Real coil = 0.0;
    Real min_source = 0.0;
};

inline AphiJoulePowerMetrics hostJoulePowerMetrics(BaseParticles &particles, const Vecd *positions,
                                                   size_t total_real_particles,
                                                   const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                                   const std::string &joule_source_name)
{
    AphiJoulePowerMetrics metrics;
    metrics.total = hostVolWeightedSum(particles, joule_source_name, total_real_particles,
                                       [](const Vecd &) { return true; });
    metrics.conductor = hostTeam7RegionJoulePower(particles, positions, total_real_particles, layout,
                                                 AphiBenchmarkMaterialRegion::Conductor, joule_source_name);
    metrics.coil = hostTeam7RegionJoulePower(particles, positions, total_real_particles, layout,
                                             AphiBenchmarkMaterialRegion::Coil, joule_source_name);

    syncVariableToHost<Real>(particles, joule_source_name);
    const Real *joule_source = particles.getVariableDataByName<Real>(joule_source_name);
    metrics.min_source = joule_source[0];
    for (size_t i = 1; i != total_real_particles; ++i)
    {
        metrics.min_source = std::min(metrics.min_source, joule_source[i]);
    }
    return metrics;
}

inline Real relativeMetricChange(Real reference, Real value)
{
    return std::abs(value - reference) / (std::abs(reference) + TinyReal);
}

struct AphiJouleLocalFieldMetrics
{
    Real l2_difference = 0.0;
    Real max_abs_difference = 0.0;
    Real conductor_l2_difference = 0.0;
    Real conductor_max_abs_difference = 0.0;
    Real conductor_band_l2_difference = 0.0;
    Real conductor_band_max_abs_difference = 0.0;
};

struct AphiJouleEmToleranceSweepEntry
{
    Real tolerance = 0.0;
    Real em_rel = 0.0;
    bool converged = false;
    UnsignedInt outer_iterations = 0;
    AphiJoulePowerMetrics joule{};
    AphiJouleLocalFieldMetrics local{};
};

struct AphiBlockJacobi8x8DiagnosticSummary
{
    size_t total_particles = 0;
    size_t fallback_count = 0;
    Real fallback_fraction = 0.0;
    Real global_min_pivot = 0.0;
    size_t conductor_fallback = 0;
    size_t coil_fallback = 0;
    size_t air_fallback = 0;
};

inline Real hostCoreMeanScalar(BaseParticles &particles, const std::string &variable_name, size_t total_real_particles,
                               Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real sum = 0.0;
    size_t count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        sum += values[i];
        ++count;
    }
    return sum / (static_cast<Real>(count) + TinyReal);
}

inline Real hostVolWeightedL2Difference(BaseParticles &particles, const std::string &variable_name,
                                        size_t total_real_particles, const Real *reference_values,
                                        const std::function<bool(const Vecd &)> &include_particle,
                                        const Vecd *positions)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_sq = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (positions != nullptr && !include_particle(positions[i]))
        {
            continue;
        }
        const Real diff = values[i] - reference_values[i];
        sum_sq += vol[i] * diff * diff;
    }
    return std::sqrt(sum_sq);
}

inline Real hostMaxAbsDifference(BaseParticles &particles, const std::string &variable_name, size_t total_real_particles,
                                 const Real *reference_values,
                                 const std::function<bool(const Vecd &)> &include_particle, const Vecd *positions)
{
    syncVariableToHost<Real>(particles, variable_name);
    const Real *values = particles.getVariableDataByName<Real>(variable_name);
    Real max_abs = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (positions != nullptr && !include_particle(positions[i]))
        {
            continue;
        }
        max_abs = std::max(max_abs, std::abs(values[i] - reference_values[i]));
    }
    return max_abs;
}

inline bool insideConductorSurfaceBand(const Vecd &position, const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                       Real band_half_width)
{
    if (!benchmark::insideBoxRegion(position, layout.conductor))
    {
        return false;
    }
    const Real dx_left = position[0] - layout.conductor.xmin;
    const Real dx_right = layout.conductor.xmax - position[0];
    const Real dy_bottom = position[1] - layout.conductor.ymin;
    const Real dy_top = layout.conductor.ymax - position[1];
    const Real dz_back = position[2] - layout.conductor.zmin;
    const Real dz_front = layout.conductor.zmax - position[2];
    const Real min_distance = std::min({dx_left, dx_right, dy_bottom, dy_top, dz_back, dz_front});
    return min_distance <= band_half_width;
}

inline AphiJouleLocalFieldMetrics hostJouleLocalFieldDifferenceMetrics(
    BaseParticles &particles, const Vecd *positions, size_t total_real_particles,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, const std::string &joule_source_name,
    const Real *reference_joule_source, Real conductor_band_half_width)
{
    const auto all_particles = [](const Vecd &) { return true; };
    const auto in_conductor = [&](const Vecd &position) {
        return team7ParticleInRegion(position, layout, AphiBenchmarkMaterialRegion::Conductor);
    };
    const auto in_conductor_band = [&](const Vecd &position) {
        return insideConductorSurfaceBand(position, layout, conductor_band_half_width);
    };

    AphiJouleLocalFieldMetrics metrics;
    metrics.l2_difference =
        hostVolWeightedL2Difference(particles, joule_source_name, total_real_particles, reference_joule_source,
                                    all_particles, positions);
    metrics.max_abs_difference = hostMaxAbsDifference(particles, joule_source_name, total_real_particles,
                                                      reference_joule_source, all_particles, positions);
    metrics.conductor_l2_difference =
        hostVolWeightedL2Difference(particles, joule_source_name, total_real_particles, reference_joule_source,
                                    in_conductor, positions);
    metrics.conductor_max_abs_difference = hostMaxAbsDifference(particles, joule_source_name, total_real_particles,
                                                                reference_joule_source, in_conductor, positions);
    metrics.conductor_band_l2_difference =
        hostVolWeightedL2Difference(particles, joule_source_name, total_real_particles, reference_joule_source,
                                      in_conductor_band, positions);
    metrics.conductor_band_max_abs_difference = hostMaxAbsDifference(particles, joule_source_name, total_real_particles,
                                                                     reference_joule_source, in_conductor_band,
                                                                     positions);
    return metrics;
}

inline AphiBlockJacobi8x8DiagnosticSummary hostBlockJacobi8x8DiagnosticSummary(
    BaseParticles &particles, const Vecd *positions, size_t total_real_particles,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, const AphiBlockJacobiDiagnosticNames &diag_names)
{
    syncVariableToHost<Real>(particles, diag_names.fallback_flag);
    syncVariableToHost<Real>(particles, diag_names.min_pivot);
    const Real *fallback_flag = particles.getVariableDataByName<Real>(diag_names.fallback_flag);
    const Real *min_pivot = particles.getVariableDataByName<Real>(diag_names.min_pivot);

    AphiBlockJacobi8x8DiagnosticSummary summary;
    summary.total_particles = total_real_particles;
    summary.global_min_pivot = min_pivot[0];
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        summary.global_min_pivot = std::min(summary.global_min_pivot, min_pivot[i]);
        if (fallback_flag[i] > Real(0.5))
        {
            ++summary.fallback_count;
            if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Conductor))
            {
                ++summary.conductor_fallback;
            }
            else if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Coil))
            {
                ++summary.coil_fallback;
            }
            else
            {
                ++summary.air_fallback;
            }
        }
    }
    summary.fallback_fraction =
        static_cast<Real>(summary.fallback_count) / (static_cast<Real>(summary.total_particles) + TinyReal);
    return summary;
}

inline std::vector<Real> hostCopyJouleSourceReference(BaseParticles &particles, const std::string &joule_source_name,
                                                      size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, joule_source_name);
    const Real *joule_source = particles.getVariableDataByName<Real>(joule_source_name);
    return std::vector<Real>(joule_source, joule_source + total_real_particles);
}

struct AphiTeam7ObservableMetrics
{
    Real em_rel = 0.0;
    UnsignedInt outer_iterations = 0;
    bool converged = false;
    Real total_joule_power = 0.0;
    Real conductor_joule_power = 0.0;
    Real coil_joule_power = 0.0;
    Real conductor_solution_block_max = 0.0;
    Real coil_solution_block_max = 0.0;
    Real conductor_E_L2 = 0.0;
    Real conductor_J_L2 = 0.0;
    Real conductor_Joule_power = 0.0;
    Real conductor_joule_max = 0.0;
};

struct AphiTeam7CenterlineProfile
{
    std::vector<Real> x_centers;
    std::vector<Real> field_norm;
};

struct AphiTeam7CanonicalCaseRunResult
{
    Real dp = 0.0;
    size_t particles = 0;
    Real core_shell = 0.0;
    benchmark::AphiTeam7LikeUnitBoxLayout layout{};
    AphiTeam7ObservableMetrics metrics{};
    AphiTeam7CenterlineProfile centerline{};
};

inline Real coreShellForTeam7PhysicalBox(Real dp, Real body_length, Real body_height, Real body_width)
{
    return std::min(Real(2.0) * dp, Real(0.125) * std::min({body_length, body_height, body_width}));
}

inline AphiTeam7CenterlineProfile hostConductorCenterlineFieldProfile(
    BaseParticles &particles, const AphiBlockNames &solution_block, const Vecd *positions, size_t total_real_particles,
    const benchmark::AphiTeam7LikeUnitBoxLayout &layout, UnsignedInt bin_count, Real yz_band_half_width)
{
    syncAphiBlockToHost(particles, solution_block);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(solution_block.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(solution_block.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(solution_block.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(solution_block.phi_imag);

    const Real y_mid = Real(0.5) * (layout.conductor.ymin + layout.conductor.ymax);
    const Real z_mid = Real(0.5) * (layout.conductor.zmin + layout.conductor.zmax);
    const Real x_min = layout.conductor.xmin;
    const Real x_max = layout.conductor.xmax;
    const Real dx = (x_max - x_min) / (static_cast<Real>(bin_count) + TinyReal);

    AphiTeam7CenterlineProfile profile;
    profile.x_centers.assign(bin_count, Real(0));
    profile.field_norm.assign(bin_count, Real(0));
    std::vector<size_t> counts(bin_count, 0);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!benchmark::insideBoxRegion(positions[i], layout.conductor))
        {
            continue;
        }
        if (std::abs(positions[i][1] - y_mid) > yz_band_half_width ||
            std::abs(positions[i][2] - z_mid) > yz_band_half_width)
        {
            continue;
        }
        const Real x = positions[i][0];
        if (x < x_min || x > x_max)
        {
            continue;
        }
        const size_t bin = static_cast<size_t>(std::min(
            static_cast<Real>(bin_count - 1), std::floor((x - x_min) / (dx + TinyReal))));
        const Real block_norm = std::sqrt(a_real[i].squaredNorm() + a_imag[i].squaredNorm() + phi_real[i] * phi_real[i] +
                                          phi_imag[i] * phi_imag[i]);
        profile.field_norm[bin] += block_norm;
        ++counts[bin];
    }

    for (UnsignedInt bin = 0; bin < bin_count; ++bin)
    {
        profile.x_centers[bin] = x_min + (static_cast<Real>(bin) + Real(0.5)) * dx;
        if (counts[bin] > 0)
        {
            profile.field_norm[bin] /= static_cast<Real>(counts[bin]);
        }
    }
    return profile;
}

inline Real hostCenterlineProfileL2RelativeDifference(const AphiTeam7CenterlineProfile &reference,
                                                      const AphiTeam7CenterlineProfile &approximate)
{
    if (reference.field_norm.size() != approximate.field_norm.size() || reference.field_norm.empty())
    {
        return Real(0);
    }
    Real sum_sq = 0.0;
    Real ref_sum_sq = 0.0;
    for (size_t i = 0; i != reference.field_norm.size(); ++i)
    {
        const Real diff = approximate.field_norm[i] - reference.field_norm[i];
        sum_sq += diff * diff;
        ref_sum_sq += reference.field_norm[i] * reference.field_norm[i];
    }
    return std::sqrt(sum_sq) / (std::sqrt(ref_sum_sq) + TinyReal);
}

inline AphiTeam7CanonicalCaseRunResult runTeam7PhysicalCanonicalCase(int ac, char *av[], Real dp_0)
{
    using benchmark::AphiTeam7CanonicalCaseSpec;
    const AphiTeam7CanonicalCaseSpec spec;

    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = coreShellForTeam7PhysicalBox(dp_0, spec.body_length, spec.body_height, spec.body_width);
    const benchmark::AphiTeam7LikeUnitBoxLayout layout =
        benchmark::buildTeam7LayoutForBox(spec.body_length, spec.body_height, spec.body_width);
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);
    const AphiJouleHeatingFieldNames joule_fields;

    AphiLhsTestBody test_body(dp_0, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);

    AphiVariableNames names;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, benchmark::AssignTeam7LikeRegionMaterialsCK> assign_material(
        test_body.body, layout, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignImpressedCurrentRhsCK> assign_coil_source(
        test_body.body, names.rhs, layout.coil, coil_current_real, coil_current_imag, spec.impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = spec.omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = spec.phi_gauge_penalty;

    AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(spec.tolerance, spec.restart_dimension, spec.max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, spec.omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);

    (void)register_joule_fields;
    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();

    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult em_result = solver.solve();

    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    AphiTeam7CanonicalCaseRunResult run_result;
    run_result.dp = dp_0;
    run_result.particles = total_real_particles;
    run_result.core_shell = core_shell;
    run_result.layout = layout;
    run_result.metrics.em_rel = em_result.final_true_relative_residual;
    run_result.metrics.outer_iterations = em_result.outer_iteration_count;
    run_result.metrics.converged = gmresConvergencePassed(em_result, spec.tolerance);
    run_result.metrics.conductor_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, spec.body_length,
                                    spec.body_height, spec.body_width, core_shell, layout.conductor);
    run_result.metrics.coil_solution_block_max = coreBlockMaxAbsNormInRegion(
        particles, names.solution, positions, total_real_particles, spec.body_length, spec.body_height, spec.body_width,
        core_shell, layout.coil);

    const AphiJoulePowerMetrics joule_metrics =
        hostJoulePowerMetrics(particles, positions, total_real_particles, layout, joule_fields.joule_heat_source);
    run_result.metrics.total_joule_power = joule_metrics.total;
    run_result.metrics.conductor_joule_power = joule_metrics.conductor;
    run_result.metrics.coil_joule_power = joule_metrics.coil;
    run_result.metrics.conductor_Joule_power = joule_metrics.conductor;

    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_real);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *joule_source = particles.getVariableDataByName<Real>(joule_fields.joule_heat_source);
    const Vecd *electric_field = particles.getVariableDataByName<Vecd>(joule_fields.electric_field_a_real);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    run_result.metrics.conductor_joule_max = 0.0;
    Real conductor_e_sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Conductor))
        {
            continue;
        }
        run_result.metrics.conductor_joule_max = std::max(run_result.metrics.conductor_joule_max, joule_source[i]);
        conductor_e_sum_squared += vol[i] * electric_field[i].squaredNorm();
    }
    run_result.metrics.conductor_E_L2 = std::sqrt(conductor_e_sum_squared);
    run_result.metrics.conductor_J_L2 = layout.conductor_material.sigma * run_result.metrics.conductor_E_L2;

    const Real yz_band = spec.centerline_yz_band_factor * dp_0;
    run_result.centerline = hostConductorCenterlineFieldProfile(particles, names.solution, positions,
                                                                total_real_particles, layout, spec.centerline_bins,
                                                                yz_band);
    return run_result;
}

enum class AphiColdCrucibleMaterialRegion
{
    Melt,
    CrucibleWall,
    Coil,
    Air
};

inline bool coldCrucibleParticleInRegion(const Vecd &position, const benchmark::AphiColdCrucibleUnitBoxLayout &layout,
                                       AphiColdCrucibleMaterialRegion region)
{
    switch (region)
    {
    case AphiColdCrucibleMaterialRegion::Melt:
        return benchmark::insideBoxRegion(position, layout.melt);
    case AphiColdCrucibleMaterialRegion::CrucibleWall:
        return benchmark::insideBoxRegion(position, layout.crucible_wall) &&
               !benchmark::insideBoxRegion(position, layout.melt);
    case AphiColdCrucibleMaterialRegion::Coil:
        return benchmark::insideAnyCoilRegion(position, layout);
    case AphiColdCrucibleMaterialRegion::Air:
    default:
        return !benchmark::insideBoxRegion(position, layout.melt) &&
               !benchmark::insideBoxRegion(position, layout.crucible_wall) &&
               !benchmark::insideAnyCoilRegion(position, layout);
    }
}

struct AphiColdCrucibleRegionCounts
{
    size_t melt = 0;
    size_t crucible_wall = 0;
    size_t coil = 0;
    size_t air = 0;
};

inline AphiColdCrucibleRegionCounts hostColdCrucibleRegionCounts(const Vecd *positions, size_t total_real_particles,
                                                               const benchmark::AphiColdCrucibleUnitBoxLayout &layout)
{
    AphiColdCrucibleRegionCounts counts;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (coldCrucibleParticleInRegion(positions[i], layout, AphiColdCrucibleMaterialRegion::Melt))
        {
            ++counts.melt;
            continue;
        }
        if (coldCrucibleParticleInRegion(positions[i], layout, AphiColdCrucibleMaterialRegion::CrucibleWall))
        {
            ++counts.crucible_wall;
            continue;
        }
        if (coldCrucibleParticleInRegion(positions[i], layout, AphiColdCrucibleMaterialRegion::Coil))
        {
            ++counts.coil;
            continue;
        }
        ++counts.air;
    }
    return counts;
}

inline Real hostColdCrucibleRegionJoulePower(BaseParticles &particles, const Vecd *positions,
                                             size_t total_real_particles,
                                             const benchmark::AphiColdCrucibleUnitBoxLayout &layout,
                                             AphiColdCrucibleMaterialRegion region, const std::string &joule_source_name)
{
    return hostVolWeightedSum(
        particles, joule_source_name, total_real_particles,
        [&](const Vecd &position) { return coldCrucibleParticleInRegion(position, layout, region); }, positions);
}

struct AphiColdCrucibleObservableMetrics
{
    Real em_rel = 0.0;
    UnsignedInt outer_iterations = 0;
    bool converged = false;
    Real total_joule_power = 0.0;
    Real melt_joule_power = 0.0;
    Real crucible_joule_power = 0.0;
    Real coil_joule_power = 0.0;
    Real melt_solution_block_max = 0.0;
    Real crucible_solution_block_max = 0.0;
    Real coil_solution_block_max = 0.0;
    Real min_joule_source = 0.0;
};

struct AphiColdCrucibleScaffoldRunResult
{
    Real dp = 0.0;
    size_t particles = 0;
    benchmark::AphiColdCrucibleUnitBoxLayout layout{};
    AphiColdCrucibleRegionCounts region_counts{};
    AphiColdCrucibleObservableMetrics metrics{};
};

inline AphiColdCrucibleScaffoldRunResult runColdCrucibleScaffoldCase(int ac, char *av[], Real dp_0)
{
    using benchmark::AphiColdCrucibleCaseSpec;
    const AphiColdCrucibleCaseSpec spec;

    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const benchmark::AphiColdCrucibleUnitBoxLayout layout =
        benchmark::buildColdCrucibleLayoutForBox(spec.body_length, spec.body_height, spec.body_width);
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);
    const AphiJouleHeatingFieldNames joule_fields;

    AphiLhsTestBody test_body(dp_0, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);

    AphiVariableNames names;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, benchmark::AssignColdCrucibleRegionMaterialsCK> assign_material(
        test_body.body, layout, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignColdCrucibleCoilSourceCK> assign_coil_source(
        test_body.body, names.rhs, layout, coil_current_real, coil_current_imag, spec.impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = spec.omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = spec.phi_gauge_penalty;

    AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(spec.tolerance, spec.restart_dimension, spec.max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, spec.omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);

    (void)register_joule_fields;
    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();

    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult em_result = solver.solve();

    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    AphiColdCrucibleScaffoldRunResult run_result;
    run_result.dp = dp_0;
    run_result.particles = total_real_particles;
    run_result.layout = layout;
    run_result.region_counts = hostColdCrucibleRegionCounts(positions, total_real_particles, layout);
    run_result.metrics.em_rel = em_result.final_true_relative_residual;
    run_result.metrics.outer_iterations = em_result.outer_iteration_count;
    run_result.metrics.converged = gmresConvergencePassed(em_result, spec.tolerance);
    run_result.metrics.melt_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, spec.body_length,
                                    spec.body_height, spec.body_width, core_shell, layout.melt);
    run_result.metrics.crucible_solution_block_max = coreBlockMaxAbsNormInRegion(
        particles, names.solution, positions, total_real_particles, spec.body_length, spec.body_height,
        spec.body_width, core_shell, layout.crucible_wall);
    run_result.metrics.coil_solution_block_max = coreBlockMaxAbsNormInRegion(
        particles, names.solution, positions, total_real_particles, spec.body_length, spec.body_height,
        spec.body_width, core_shell, layout.coil_left);

    run_result.metrics.total_joule_power =
        hostVolWeightedSum(particles, joule_fields.joule_heat_source, total_real_particles,
                           [](const Vecd &) { return true; });
    run_result.metrics.melt_joule_power = hostColdCrucibleRegionJoulePower(
        particles, positions, total_real_particles, layout, AphiColdCrucibleMaterialRegion::Melt,
        joule_fields.joule_heat_source);
    run_result.metrics.crucible_joule_power = hostColdCrucibleRegionJoulePower(
        particles, positions, total_real_particles, layout, AphiColdCrucibleMaterialRegion::CrucibleWall,
        joule_fields.joule_heat_source);
    run_result.metrics.coil_joule_power = hostColdCrucibleRegionJoulePower(
        particles, positions, total_real_particles, layout, AphiColdCrucibleMaterialRegion::Coil,
        joule_fields.joule_heat_source);

    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    const Real *joule_source = particles.getVariableDataByName<Real>(joule_fields.joule_heat_source);
    run_result.metrics.min_joule_source = joule_source[0];
    for (size_t i = 1; i != total_real_particles; ++i)
    {
        run_result.metrics.min_joule_source = std::min(run_result.metrics.min_joule_source, joule_source[i]);
    }
    return run_result;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_BENCHMARK_HELPERS_H
