#ifndef APHI_DIVERGENCE_FREE_A_MMS_HELPERS_H
#define APHI_DIVERGENCE_FREE_A_MMS_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_assemble_lhs_debug_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_gmres_robustness_sweep_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_eta_sweep_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_inner_a_divergence_penalty_gate_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

enum class AphiDivergenceFreeAFieldKind
{
    /** A=(-y/2,x/2,0), analytic divA=0; matches curl_a_manufactured. */
    Linear2D,
    /** A=curl(0,0,psi), psi=sin(pi x)sin(pi y)sin(pi z). */
    Sinusoidal3DCurlPsi
};

inline Vecd divergenceFreeARealField(AphiDivergenceFreeAFieldKind kind, Real x, Real y, Real z)
{
    if (kind == AphiDivergenceFreeAFieldKind::Linear2D)
    {
        return Vecd(-0.5 * y, 0.5 * x, 0.0);
    }
    const Real pi = Pi;
    const Real sx = std::sin(pi * x);
    const Real cx = std::cos(pi * x);
    const Real sy = std::sin(pi * y);
    const Real cy = std::cos(pi * y);
    const Real sz = std::sin(pi * z);
    return Vecd(pi * sx * cy * sz, -pi * cx * sy * sz, 0.0);
}

class AssignDivergenceFreeAphiFieldCK : public LocalDynamics
{
  public:
    AssignDivergenceFreeAphiFieldCK(SPHBody &sph_body, const AphiBlockNames &block_names,
                                      AphiDivergenceFreeAFieldKind field_kind, bool use_nonzero_phi = false)
        : LocalDynamics(sph_body), field_kind_(field_kind), use_nonzero_phi_(use_nonzero_phi),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : field_kind_(encloser.field_kind_), use_nonzero_phi_(encloser.use_nonzero_phi_),
              position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            a_real_[index_i] = divergenceFreeARealField(field_kind_, position[0], position[1], position[2]);
            a_imag_[index_i] = Vecd(0.0, 0.0, 0.0);
            if (use_nonzero_phi_)
            {
                const Real pi = Pi;
                const Real x = position[0];
                const Real y = position[1];
                const Real z = position[2];
                phi_real_[index_i] = 0.1 * std::sin(pi * x) * std::cos(pi * y);
                phi_imag_[index_i] = 0.05 * std::cos(pi * x) * std::sin(pi * z);
            }
            else
            {
                phi_real_[index_i] = 0.0;
                phi_imag_[index_i] = 0.0;
            }
        }

      protected:
        AphiDivergenceFreeAFieldKind field_kind_;
        bool use_nonzero_phi_;
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    AphiDivergenceFreeAFieldKind field_kind_;
    bool use_nonzero_phi_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

struct AphiDivergenceFreeADivADiagnosticMetrics
{
    AphiDivergenceFreeAFieldKind field_kind = AphiDivergenceFreeAFieldKind::Linear2D;
    AphiDivAReductionMetrics b_corrected{};
    AphiDivAReductionMetrics pairwise{};
    size_t core_particles = 0;
    size_t total_particles = 0;
};

inline void registerDivergenceFreeExactBlock(BaseParticles &particles, const AphiBlockNames &exact_block)
{
    particles.registerStateVariable<Vecd>(exact_block.a_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(exact_block.a_imag, ZeroData<Vecd>::value);
    particles.registerStateVariable<Real>(exact_block.phi_real, Real(0));
    particles.registerStateVariable<Real>(exact_block.phi_imag, Real(0));
}

inline AphiDivergenceFreeADivADiagnosticMetrics runDivergenceFreeADivADiagnostic(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, AphiDivergenceFreeAFieldKind field_kind)
{
    AphiDivergenceFreeADivADiagnosticMetrics metrics;
    metrics.field_kind = field_kind;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivergenceFreeAphiFieldCK> assign_field(
        test_body.body, names.solution, field_kind);

    initialize_aphi.exec();
    set_material.exec();
    assign_field.exec();
    test_body.updateRelations();

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names, AphiDivADiagnosticMode::BCorrectedTrace);
    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();
    metrics.b_corrected = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, {},
                                                       AphiDivADiagnosticMode::BCorrectedTrace);

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names, AphiDivADiagnosticMode::PairwiseUncorrected);
    metrics.pairwise = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, {},
                                                    AphiDivADiagnosticMode::PairwiseUncorrected);
    metrics.core_particles = metrics.pairwise.particle_count;
    metrics.total_particles = total_real_particles;
    return metrics;
}

struct AphiDivergenceFreeAMMSRow
{
    AphiDivergenceFreeAFieldKind field_kind = AphiDivergenceFreeAFieldKind::Sinusoidal3DCurlPsi;
    Real eta_a = 0.0;
    Real lambda_a = 0.0;
    bool consistency_rhs = true;
    bool converged = false;
    Real rhs_norm = 0.0;
    Real true_rel = 0.0;
    Real a_linf_error = 0.0;
    Real phi_linf_error = 0.0;
    Real block_linf_error = 0.0;
    AphiDivAReductionMetrics div_a{};
    Real joule_power = 0.0;
    Real E_L2 = 0.0;
    Real J_L2 = 0.0;
};

inline AphiDivergenceFreeAMMSRow runManufacturedSeparableAphiMMSRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, bool consistency_rhs,
    UnsignedInt restart_dimension = 50, UnsignedInt max_outer_iterations = 100)
{
    AphiDivergenceFreeAMMSRow row;
    row.eta_a = eta_a;
    row.consistency_rhs = consistency_rhs;
    row.lambda_a = eta_a > TinyReal ? lambdaAFromEtaA(eta_a, scale_metrics) : 0.0;

    const AphiBlockNames exact_block{"ExactAReal", "ExactAImag", "ExactPhiReal", "ExactPhiImag"};
    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;
    const AphiJouleHeatingFieldNames joule_fields;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    registerDivergenceFreeExactBlock(test_body.body.getBaseParticles(), exact_block);

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiManufacturedFieldsCK> assign_exact(
        test_body.body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact(test_body.body, exact_block, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions rhs_options;
    rhs_options.omega = omega;
    rhs_options.use_phi_gauge_penalty = true;
    rhs_options.phi_gauge_penalty = phi_gauge_penalty;
    rhs_options.use_a_divergence_penalty = consistency_rhs ? (eta_a > TinyReal) : false;
    rhs_options.a_divergence_penalty = consistency_rhs ? row.lambda_a : 0.0;

    AphiLhsAssemblyOptions solve_options = rhs_options;
    solve_options.use_a_divergence_penalty = eta_a > TinyReal;
    solve_options.a_divergence_penalty = row.lambda_a;

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_exact(test_body.body, test_body.inner(), names.solution, names.lhs,
                                                             names.material, omega, rhs_options);
    const AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, solve_options,
                                                        solver_options);

    initialize_aphi.exec();
    set_material.exec();
    assign_exact.exec();
    test_body.updateRelations();
    copy_exact.exec();
    apply_exact.exec();
    copy_lhs_to_rhs.exec();
    zero_solution.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncAphiBlockToHost(particles, names.rhs);
    row.rhs_norm = hostBlockNorm(particles, names.rhs, total_real_particles);

    const AphiMatrixFreeSolverResult result = solver.solve();
    row.converged = result.converged;
    row.true_rel = result.final_true_relative_residual;

    (void)register_joule_fields;
    execInnerJoulePostProcess(test_body, names, omega, joule_fields);
    syncAphiBlockToHost(particles, names.solution);
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    row.block_linf_error = coreBlockMaxAbsDifference(particles, names.solution, exact_block, positions,
                                                     total_real_particles, body_length, body_height, body_width,
                                                     core_shell);
    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names);
    test_body.updateRelations();
    row.div_a = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles);
    const auto in_core = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    const AphiRegionalElectromagneticObservables observables = hostRegionalElectromagneticObservablesFromJouleFields(
        particles, names, joule_fields, positions, total_real_particles, body_length, body_height, body_width,
        core_shell, in_core, sigma);
    row.joule_power = observables.Joule_power;
    row.E_L2 = observables.E_L2;
    row.J_L2 = observables.J_L2;
    return row;
}

inline AphiDivergenceFreeAMMSRow runDivergenceFreeAMMSRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real tolerance, AphiDivergenceFreeAFieldKind field_kind,
    bool consistency_rhs, UnsignedInt restart_dimension = 50, UnsignedInt max_outer_iterations = 100)
{
    AphiDivergenceFreeAMMSRow row;
    row.field_kind = field_kind;
    row.eta_a = eta_a;
    row.consistency_rhs = consistency_rhs;
    row.lambda_a = eta_a > TinyReal ? lambdaAFromEtaA(eta_a, scale_metrics) : 0.0;

    const AphiBlockNames exact_block{"ExactAReal", "ExactAImag", "ExactPhiReal", "ExactPhiImag"};

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    registerDivergenceFreeExactBlock(test_body.body.getBaseParticles(), exact_block);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivergenceFreeAphiFieldCK> assign_exact(
        test_body.body, names.solution, field_kind, true);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact(test_body.body, exact_block, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions rhs_options;
    rhs_options.omega = omega;
    rhs_options.use_phi_gauge_penalty = true;
    rhs_options.phi_gauge_penalty = phi_gauge_penalty;
    rhs_options.use_a_divergence_penalty = consistency_rhs ? (eta_a > TinyReal) : false;
    rhs_options.a_divergence_penalty = consistency_rhs ? row.lambda_a : 0.0;

    AphiLhsAssemblyOptions solve_options = rhs_options;
    solve_options.use_a_divergence_penalty = eta_a > TinyReal;
    solve_options.a_divergence_penalty = row.lambda_a;

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble_exact(test_body.body, test_body.inner(), names,
                                                                           rhs_options);
    const AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, solve_options,
                                                        solver_options);

    initialize_aphi.exec();
    set_material.exec();
    assign_exact.exec();
    test_body.updateRelations();
    copy_exact.exec();
    assemble_exact.exec();
    copy_lhs_to_rhs.exec();
    zero_solution.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncAphiBlockToHost(particles, names.rhs);
    syncAphiBlockToHost(particles, names.lhs);
    row.rhs_norm = hostBlockNorm(particles, names.rhs, total_real_particles);

    const AphiMatrixFreeSolverResult result = solver.solve();
    row.converged = result.converged;
    row.true_rel = result.final_true_relative_residual;

    syncAphiBlockToHost(particles, names.solution);
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    row.block_linf_error = coreBlockMaxAbsDifference(particles, names.solution, exact_block, positions,
                                                     total_real_particles, body_length, body_height, body_width,
                                                     core_shell);
    syncAphiBlockToHost(particles, names.solution);
    syncAphiBlockToHost(particles, exact_block);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(names.solution.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(names.solution.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(names.solution.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(names.solution.phi_imag);
    const Vecd *exact_a_real = particles.getVariableDataByName<Vecd>(exact_block.a_real);
    const Vecd *exact_a_imag = particles.getVariableDataByName<Vecd>(exact_block.a_imag);
    const Real *exact_phi_real = particles.getVariableDataByName<Real>(exact_block.phi_real);
    const Real *exact_phi_imag = particles.getVariableDataByName<Real>(exact_block.phi_imag);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        row.a_linf_error = std::max(row.a_linf_error, (a_real[i] - exact_a_real[i]).norm());
        row.a_linf_error = std::max(row.a_linf_error, (a_imag[i] - exact_a_imag[i]).norm());
        row.phi_linf_error = std::max(row.phi_linf_error, std::abs(phi_real[i] - exact_phi_real[i]));
        row.phi_linf_error = std::max(row.phi_linf_error, std::abs(phi_imag[i] - exact_phi_imag[i]));
    }

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names);
    test_body.updateRelations();
    row.div_a = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles);
    return row;
}

inline bool divergenceFreeAMMSInvarianceRowPassed(const AphiDivergenceFreeAMMSRow &row, Real tolerance,
                                                  Real min_rhs_norm = 1.0e-8)
{
    if (row.rhs_norm < min_rhs_norm)
    {
        return false;
    }
    return row.converged && std::isfinite(row.true_rel) && row.true_rel <= tolerance;
}

inline bool divergenceFreeAMMSRowPassed(const AphiDivergenceFreeAMMSRow &row, Real tolerance, Real max_block_linf_error,
                                        Real min_rhs_norm = 1.0e-8)
{
    if (row.rhs_norm < min_rhs_norm)
    {
        return false;
    }
    if (!row.converged || !std::isfinite(row.true_rel) || row.true_rel > tolerance)
    {
        return false;
    }
    return row.block_linf_error <= max_block_linf_error;
}

inline void printDivergenceFreeADivADiagnosticMetrics(const char *test_name,
                                                      const AphiDivergenceFreeADivADiagnosticMetrics &metrics)
{
    const char *field_name =
        metrics.field_kind == AphiDivergenceFreeAFieldKind::Linear2D ? "Linear2D" : "Sinusoidal3DCurlPsi";
    std::cout << test_name << " field=" << field_name << " core_particles=" << metrics.core_particles
              << " divA_B_L2=" << metrics.b_corrected.div_a_L2 << " divA_B_rel=" << metrics.b_corrected.div_a_relative
              << " divA_pairwise_L2=" << metrics.pairwise.div_a_L2
              << " divA_pairwise_rel=" << metrics.pairwise.div_a_relative << std::endl;
}

struct AphiDivergenceFreeADpRefinementRow
{
    Real dp = 0.0;
    size_t total_particles = 0;
    Real linear2d_pairwise_div_a_rel = 0.0;
    Real linear2d_b_div_a_rel = 0.0;
    Real mms_eta_a = 0.0;
    Real mms_true_rel = 0.0;
    Real mms_block_linf_error = 0.0;
    Real mms_div_a_rel = 0.0;
};

inline AphiDivergenceFreeADpRefinementRow runDivergenceFreeADpRefinementRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, Real omega, Real phi_gauge_penalty, Real eta_a, const AphiCoreOperatorScaleMetrics &scale_metrics,
    Real tolerance)
{
    const Real boundary_width = 3.0 * dp_0;
    AphiDivergenceFreeADpRefinementRow row;
    row.dp = dp_0;
    row.mms_eta_a = eta_a;

    const AphiDivergenceFreeADivADiagnosticMetrics linear2d = runDivergenceFreeADivADiagnostic(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu,
        AphiDivergenceFreeAFieldKind::Linear2D);
    row.total_particles = linear2d.total_particles;
    row.linear2d_pairwise_div_a_rel = linear2d.pairwise.div_a_relative;
    row.linear2d_b_div_a_rel = linear2d.b_corrected.div_a_relative;

    const AphiDivergenceFreeAMMSRow mms_row = runManufacturedSeparableAphiMMSRow(
        ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, sigma, nu, omega,
        phi_gauge_penalty, eta_a, scale_metrics, tolerance, true);
    row.mms_true_rel = mms_row.true_rel;
    row.mms_block_linf_error = mms_row.block_linf_error;
    row.mms_div_a_rel = mms_row.div_a.div_a_relative;
    return row;
}

inline void printDivergenceFreeADpRefinementRow(const char *test_name, const AphiDivergenceFreeADpRefinementRow &row)
{
    std::cout << test_name << " dp=" << row.dp << " particles=" << row.total_particles << " eta_A=" << row.mms_eta_a
              << " linear2d_pairwise_divA_rel=" << row.linear2d_pairwise_div_a_rel
              << " linear2d_B_divA_rel=" << row.linear2d_b_div_a_rel << " mms_true_rel=" << row.mms_true_rel
              << " mms_block_Linf_err=" << row.mms_block_linf_error << " mms_divA_rel=" << row.mms_div_a_rel
              << std::endl;
}

inline bool divergenceFreeADpRefinementPassed(const std::vector<AphiDivergenceFreeADpRefinementRow> &rows,
                                              Real tolerance, Real max_linear2d_pairwise_div_a_rel,
                                              Real max_mms_block_linf_per_dp)
{
    if (rows.empty())
    {
        return false;
    }
    for (const AphiDivergenceFreeADpRefinementRow &row : rows)
    {
        if (row.linear2d_pairwise_div_a_rel > max_linear2d_pairwise_div_a_rel)
        {
            return false;
        }
        if (!std::isfinite(row.mms_true_rel) || row.mms_true_rel > tolerance)
        {
            return false;
        }
        if (row.mms_block_linf_error > max_mms_block_linf_per_dp)
        {
            return false;
        }
    }
    return true;
}

inline void printDivergenceFreeAMMSRow(const char *test_name, const AphiDivergenceFreeAMMSRow &row)
{
    const char *field_name =
        row.field_kind == AphiDivergenceFreeAFieldKind::Linear2D ? "Linear2D" : "ManufacturedSeparable";
    std::cout << test_name << " field=" << field_name << " mode=" << (row.consistency_rhs ? "consistency" : "invariance")
              << " eta_A=" << row.eta_a << " lambda_A=" << row.lambda_a << " rhs_norm=" << row.rhs_norm
              << " converged=" << (row.converged ? 1 : 0)
              << " true_rel=" << row.true_rel << " A_Linf_err=" << row.a_linf_error
              << " phi_Linf_err=" << row.phi_linf_error << " block_Linf_err=" << row.block_linf_error
              << " div_A_rel=" << row.div_a.div_a_relative << " div_A_L2=" << row.div_a.div_a_L2 << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_DIVERGENCE_FREE_A_MMS_HELPERS_H
