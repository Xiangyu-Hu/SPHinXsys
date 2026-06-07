#ifndef APHI_DIVERGENCE_FREE_NONZERO_RHS_MMS_HELPERS_H
#define APHI_DIVERGENCE_FREE_NONZERO_RHS_MMS_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_divergence_free_a_mms_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_a_gauge_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_curl_b_dual_track_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Non-degenerate divergence-free validation fields (Stage 10.10). */
enum class AphiDivFreeValidationFieldKind
{
    Linear2D,
    /** A=(0,0,sin(pi x)sin(pi y)), divA=0, laplace Az != 0. */
    Az2D,
    /** Cross-sine 3D field, divA=0. */
    CrossSine3D,
    Sinusoidal3DCurlPsi
};

inline const char *divFreeValidationFieldName(AphiDivFreeValidationFieldKind kind)
{
    switch (kind)
    {
    case AphiDivFreeValidationFieldKind::Linear2D:
        return "Linear2D";
    case AphiDivFreeValidationFieldKind::Az2D:
        return "Az2D";
    case AphiDivFreeValidationFieldKind::CrossSine3D:
        return "CrossSine3D";
    default:
        return "Sinusoidal3DCurlPsi";
    }
}

inline Vecd divFreeValidationARealField(AphiDivFreeValidationFieldKind kind, Real x, Real y, Real z)
{
    if (kind == AphiDivFreeValidationFieldKind::Linear2D)
    {
        return divergenceFreeARealField(AphiDivergenceFreeAFieldKind::Linear2D, x, y, z);
    }
    if (kind == AphiDivFreeValidationFieldKind::Sinusoidal3DCurlPsi)
    {
        return divergenceFreeARealField(AphiDivergenceFreeAFieldKind::Sinusoidal3DCurlPsi, x, y, z);
    }
    const Real pi = Pi;
    const Real sx = std::sin(pi * x);
    const Real sy = std::sin(pi * y);
    const Real sz = std::sin(pi * z);
    if (kind == AphiDivFreeValidationFieldKind::Az2D)
    {
        return Vecd(0.0, 0.0, sx * sy);
    }
    return Vecd(sy * sz, sz * sx, sx * sy);
}

inline Vecd divFreeValidationBRealField(AphiDivFreeValidationFieldKind kind, Real x, Real y, Real z)
{
    const Real pi = Pi;
    const Real sx = std::sin(pi * x);
    const Real sy = std::sin(pi * y);
    const Real sz = std::sin(pi * z);
    const Real cx = std::cos(pi * x);
    const Real cy = std::cos(pi * y);
    const Real cz = std::cos(pi * z);
    if (kind == AphiDivFreeValidationFieldKind::Linear2D)
    {
        return Vecd(0.0, 0.0, 1.0);
    }
    if (kind == AphiDivFreeValidationFieldKind::Az2D)
    {
        return Vecd(pi * sx * cy, -pi * cx * sy, 0.0);
    }
    if (kind == AphiDivFreeValidationFieldKind::CrossSine3D)
    {
        return Vecd(pi * sx * cy - pi * cz * sx, pi * sy * cz - pi * cx * sy, pi * sz * cx - pi * cy * sz);
    }
    return Vecd(0.0, 0.0, 0.0);
}

inline void execBodyCurlBFromADiagnostic(SPHBody &body, Inner<> &inner, const AphiVariableNames &names,
                                         const std::string &b_real_name, const std::string &b_imag_name)
{
    execBodyCurlBFromADiagnostic(body, inner, names, b_real_name, b_imag_name, AphiBCurlDiagnosticMode::BCorrectedGrad);
}

inline Real hostCoreDivFreeBRelativeError(
    BaseParticles &particles, const std::string &b_real_name, const std::string &b_imag_name, const Vecd *positions,
    size_t total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell,
    AphiDivFreeValidationFieldKind field_kind)
{
    syncVariableToHost<Vecd>(particles, b_real_name);
    syncVariableToHost<Vecd>(particles, b_imag_name);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *b_real = particles.getVariableDataByName<Vecd>(b_real_name);
    const Vecd *b_imag = particles.getVariableDataByName<Vecd>(b_imag_name);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real error_squared = 0.0;
    Real reference_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        const Vecd exact_b = divFreeValidationBRealField(field_kind, positions[i][0], positions[i][1], positions[i][2]);
        const Vecd err_re = b_real[i] - exact_b;
        const Vecd err_im = b_imag[i];
        error_squared += vol[i] * (err_re.squaredNorm() + err_im.squaredNorm());
        reference_squared += vol[i] * exact_b.squaredNorm();
    }
    return std::sqrt(error_squared) / (std::sqrt(reference_squared) + TinyReal);
}

class AssignDivFreeValidationAphiFieldCK : public LocalDynamics
{
  public:
    AssignDivFreeValidationAphiFieldCK(SPHBody &sph_body, const AphiBlockNames &block_names,
                                       AphiDivFreeValidationFieldKind field_kind)
        : LocalDynamics(sph_body), field_kind_(field_kind),
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
            : field_kind_(encloser.field_kind_),
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
            a_real_[index_i] =
                divFreeValidationARealField(field_kind_, position[0], position[1], position[2]);
            a_imag_[index_i] = Vecd(0.0, 0.0, 0.0);
            phi_real_[index_i] = 0.0;
            phi_imag_[index_i] = 0.0;
        }

      protected:
        AphiDivFreeValidationFieldKind field_kind_;
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    AphiDivFreeValidationFieldKind field_kind_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

struct AphiDivFreeValidationDivADiagnosticMetrics
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real dp = 0.0;
    Real core_shell = 0.0;
    AphiDivAReductionMetrics b_corrected{};
    AphiDivAReductionMetrics pairwise{};
    size_t core_particles = 0;
    size_t total_particles = 0;
};

inline AphiDivFreeValidationDivADiagnosticMetrics runDivFreeValidationDivADiagnostic(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, AphiDivFreeValidationFieldKind field_kind)
{
    const Real boundary_width = 3.0 * dp_0;
    AphiDivFreeValidationDivADiagnosticMetrics metrics;
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

struct AphiNonDegenerateDivFreeMMSRow
{
    AphiDivFreeValidationFieldKind field_kind = AphiDivFreeValidationFieldKind::Az2D;
    Real eta_a = 0.0;
    Real lambda_a = 0.0;
    bool consistency_rhs = true;
    bool converged = false;
    Real rhs_norm = 0.0;
    Real true_rel = 0.0;
    Real block_linf_error = 0.0;
    AphiDivAReductionMetrics div_a_pairwise{};
    AphiDivAReductionMetrics div_a_b_corrected{};
    Real joule_power = 0.0;
    Real E_L2 = 0.0;
    Real J_L2 = 0.0;
    Real exact_joule_power = 0.0;
    Real exact_E_L2 = 0.0;
    Real exact_J_L2 = 0.0;
    Real joule_error_vs_exact = 0.0;
    Real E_L2_error_vs_exact = 0.0;
    Real J_L2_error_vs_exact = 0.0;
    Real E_combined_L2 = 0.0;
    Real exact_E_combined_L2 = 0.0;
    Real E_combined_error_vs_exact = 0.0;
    Real J_combined_L2 = 0.0;
    Real exact_J_combined_L2 = 0.0;
    Real J_combined_error_vs_exact = 0.0;
    Real b_error_vs_exact = 0.0;
    Real h_error_vs_exact = 0.0;
    Real div_a_core_grad_den_rel = 0.0;
    Real div_a_boundary_grad_den_rel = 0.0;
};

inline Real hostCoreCombinedCurrentDensityL2(
    BaseParticles &particles, const AphiJouleHeatingFieldNames &joule_fields, const Vecd *positions,
    size_t total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Vecd>(particles, joule_fields.current_density_real);
    syncVariableToHost<Vecd>(particles, joule_fields.current_density_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *j_real = particles.getVariableDataByName<Vecd>(joule_fields.current_density_real);
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(joule_fields.current_density_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        sum_squared += vol[i] * (j_real[i].squaredNorm() + j_imag[i].squaredNorm());
    }
    return std::sqrt(sum_squared);
}

inline Real hostCoreCombinedElectricFieldL2(
    BaseParticles &particles, const AphiJouleHeatingFieldNames &joule_fields, const Vecd *positions,
    size_t total_real_particles, Real body_length, Real body_height, Real body_width, Real core_shell)
{
    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_real);
    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *e_real = particles.getVariableDataByName<Vecd>(joule_fields.electric_field_a_real);
    const Vecd *e_imag = particles.getVariableDataByName<Vecd>(joule_fields.electric_field_a_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_squared = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        sum_squared += vol[i] * (e_real[i].squaredNorm() + e_imag[i].squaredNorm());
    }
    return std::sqrt(sum_squared);
}

inline AphiRegionalElectromagneticObservables hostCoreObservablesFromJouleFields(
    BaseParticles &particles, const AphiVariableNames &names, const AphiJouleHeatingFieldNames &joule_fields,
    const Vecd *positions, size_t total_real_particles, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real sigma)
{
    const auto in_core = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    return hostRegionalElectromagneticObservablesFromJouleFields(particles, names, joule_fields, positions,
                                                                 total_real_particles, body_length, body_height,
                                                                 body_width, core_shell, in_core, sigma);
}

inline AphiNonDegenerateDivFreeMMSRow runNonDegenerateDivFreeMMSRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real core_shell, Real sigma,
    Real nu, Real omega, Real phi_gauge_penalty, Real eta_a, const AphiCoreOperatorScaleMetrics &scale_metrics,
    Real tolerance, AphiDivFreeValidationFieldKind field_kind, bool consistency_rhs,
    UnsignedInt restart_dimension = 50, UnsignedInt max_outer_iterations = 100)
{
    AphiNonDegenerateDivFreeMMSRow row;
    row.field_kind = field_kind;
    row.eta_a = eta_a;
    row.consistency_rhs = consistency_rhs;
    row.lambda_a = eta_a > TinyReal ? lambdaAFromEtaA(eta_a, scale_metrics) : 0.0;

    const Real boundary_width = 3.0 * dp_0;
    const AphiBlockNames exact_block{"ExactAReal", "ExactAImag", "ExactPhiReal", "ExactPhiImag"};

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);
    AphiVariableNames names;
    const AphiJouleHeatingFieldNames joule_fields;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    registerDivergenceFreeExactBlock(test_body.body.getBaseParticles(), exact_block);
    const std::string b_real_name = "DiagnosticBReal";
    const std::string b_imag_name = "DiagnosticBImag";
    test_body.body.getBaseParticles().registerStateVariable<Vecd>(b_real_name, ZeroData<Vecd>::value);
    test_body.body.getBaseParticles().registerStateVariable<Vecd>(b_imag_name, ZeroData<Vecd>::value);

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivFreeValidationAphiFieldCK> assign_exact(
        test_body.body, names.solution, field_kind);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact(test_body.body, exact_block, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> backup_solution(test_body.body, names.t, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> restore_solution(test_body.body, names.solution, names.t);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact_to_solution(test_body.body, names.solution, exact_block);
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

    (void)register_joule_fields;
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

    syncAphiBlockToHost(particles, names.solution);
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    row.block_linf_error = coreBlockMaxAbsDifference(particles, names.solution, exact_block, positions,
                                                     total_real_particles, body_length, body_height, body_width,
                                                     core_shell);

    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names, AphiDivADiagnosticMode::PairwiseUncorrected);
    row.div_a_pairwise = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, {},
                                                      AphiDivADiagnosticMode::PairwiseUncorrected);
    execBodyDivADiagnosticPipeline(test_body.body, test_body.inner(), names, AphiDivADiagnosticMode::BCorrectedTrace);
    row.div_a_b_corrected = hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, {},
                                                         AphiDivADiagnosticMode::BCorrectedTrace);
    row.div_a_core_grad_den_rel =
        hostCoreDivARelative(particles, names, total_real_particles, body_length, body_height, body_width, core_shell);
    const auto in_boundary = [&](const Vecd &position) {
        return !isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    row.div_a_boundary_grad_den_rel =
        hostBodyDivAReductionMetrics(particles, names, positions, total_real_particles, in_boundary,
                                     AphiDivADiagnosticMode::BCorrectedTrace)
            .div_a_relative;

    execBodyCurlBFromADiagnostic(test_body.body, test_body.inner(), names, b_real_name, b_imag_name);
    row.b_error_vs_exact = hostCoreDivFreeBRelativeError(particles, b_real_name, b_imag_name, positions,
                                                         total_real_particles, body_length, body_height, body_width,
                                                         core_shell, field_kind);
    row.h_error_vs_exact = row.b_error_vs_exact;

    execInnerJoulePostProcess(test_body, names, omega, joule_fields);
    const AphiRegionalElectromagneticObservables solved_observables =
        hostCoreObservablesFromJouleFields(particles, names, joule_fields, positions, total_real_particles,
                                           body_length, body_height, body_width, core_shell, sigma);
    row.joule_power = solved_observables.Joule_power;
    row.E_L2 = solved_observables.E_L2;
    row.J_L2 = solved_observables.J_L2;
    row.E_combined_L2 = hostCoreCombinedElectricFieldL2(particles, joule_fields, positions, total_real_particles,
                                                        body_length, body_height, body_width, core_shell);
    row.J_combined_L2 = hostCoreCombinedCurrentDensityL2(particles, joule_fields, positions, total_real_particles,
                                                       body_length, body_height, body_width, core_shell);

    backup_solution.exec();
    copy_exact_to_solution.exec();
    test_body.updateRelations();
    execInnerJoulePostProcess(test_body, names, omega, joule_fields);
    const AphiRegionalElectromagneticObservables exact_observables =
        hostCoreObservablesFromJouleFields(particles, names, joule_fields, positions, total_real_particles,
                                           body_length, body_height, body_width, core_shell, sigma);
    restore_solution.exec();
    test_body.updateRelations();

    row.exact_joule_power = exact_observables.Joule_power;
    row.exact_E_L2 = exact_observables.E_L2;
    row.exact_J_L2 = exact_observables.J_L2;
    row.exact_E_combined_L2 = hostCoreCombinedElectricFieldL2(particles, joule_fields, positions, total_real_particles,
                                                              body_length, body_height, body_width, core_shell);
    row.exact_J_combined_L2 = hostCoreCombinedCurrentDensityL2(particles, joule_fields, positions, total_real_particles,
                                                               body_length, body_height, body_width, core_shell);
    row.joule_error_vs_exact = relativeMetricChange(row.exact_joule_power, row.joule_power);
    row.E_L2_error_vs_exact = relativeMetricChange(row.exact_E_L2, row.E_L2);
    row.J_L2_error_vs_exact = relativeMetricChange(row.exact_J_L2, row.J_L2);
    row.E_combined_error_vs_exact = relativeMetricChange(row.exact_E_combined_L2, row.E_combined_L2);
    row.J_combined_error_vs_exact = relativeMetricChange(row.exact_J_combined_L2, row.J_combined_L2);
    return row;
}

inline void printDivFreeValidationDivADiagnosticMetrics(const char *test_name,
                                                      const AphiDivFreeValidationDivADiagnosticMetrics &metrics)
{
    std::cout << test_name << " field=" << divFreeValidationFieldName(metrics.field_kind) << " dp=" << metrics.dp
              << " core_shell=" << metrics.core_shell << " core_particles=" << metrics.core_particles
              << " pairwise_divA_L2=" << metrics.pairwise.div_a_L2
              << " pairwise_divA_rel=" << metrics.pairwise.div_a_relative
              << " pairwise_gradA_L2=" << metrics.pairwise.grad_a_L2
              << " B_divA_L2=" << metrics.b_corrected.div_a_L2
              << " B_divA_rel=" << metrics.b_corrected.div_a_relative
              << " B_gradA_L2=" << metrics.b_corrected.grad_a_L2 << std::endl;
}

inline void printNonDegenerateDivFreeMMSRow(const char *test_name, const AphiNonDegenerateDivFreeMMSRow &row)
{
    std::cout << test_name << " field=" << divFreeValidationFieldName(row.field_kind)
              << " mode=" << (row.consistency_rhs ? "consistency" : "invariance") << " eta_A=" << row.eta_a
              << " lambda_A=" << row.lambda_a << " rhs_norm=" << row.rhs_norm << " converged=" << (row.converged ? 1 : 0)
              << " true_rel=" << row.true_rel << " block_Linf_err=" << row.block_linf_error
              << " divA_pairwise_rel=" << row.div_a_pairwise.div_a_relative
              << " divA_B_rel=" << row.div_a_b_corrected.div_a_relative << " joule_power=" << row.joule_power
              << " exact_joule_power=" << row.exact_joule_power               << " joule_err_vs_exact=" << row.joule_error_vs_exact
              << " E_L2_err_vs_exact=" << row.E_L2_error_vs_exact
              << " E_combined_err_vs_exact=" << row.E_combined_error_vs_exact
              << " J_combined_err_vs_exact=" << row.J_combined_error_vs_exact
              << " B_err_vs_exact=" << row.b_error_vs_exact << " H_err_vs_exact=" << row.h_error_vs_exact
              << " divA_core_gradDen_rel=" << row.div_a_core_grad_den_rel
              << " divA_boundary_gradDen_rel=" << row.div_a_boundary_grad_den_rel
              << " J_L2_err_vs_exact=" << row.J_L2_error_vs_exact << std::endl;
}

inline void printInnerEtaObservableMatrixSummary(const char *test_name,
                                                 const StdVec<AphiNonDegenerateDivFreeMMSRow> &rows)
{
    std::cout << test_name << " eta_observable_matrix field eta_A mode conv true_rel A_err B_err E_comb_err J_comb_err "
                 "Joule_err core_divA boundary_divA"
              << std::endl;
    for (const AphiNonDegenerateDivFreeMMSRow &row : rows)
    {
        std::cout << test_name << " eta_observable_matrix " << divFreeValidationFieldName(row.field_kind) << " "
                  << row.eta_a << " " << (row.consistency_rhs ? "consistency" : "invariance") << " "
                  << (row.converged ? 1 : 0) << " " << row.true_rel << " " << row.block_linf_error << " "
                  << row.b_error_vs_exact << " " << row.E_combined_error_vs_exact << " "
                  << row.J_combined_error_vs_exact << " " << row.joule_error_vs_exact << " "
                  << row.div_a_core_grad_den_rel << " " << row.div_a_boundary_grad_den_rel << std::endl;
    }
}

inline bool nonDegenerateDivFreeMMSConsistencyRowPassed(const AphiNonDegenerateDivFreeMMSRow &row, Real tolerance,
                                                        Real max_block_linf_error, Real min_rhs_norm = 1.0e-6)
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

inline bool nonDegenerateDivFreeMMSInvarianceRowPassed(const AphiNonDegenerateDivFreeMMSRow &row, Real tolerance,
                                                       Real max_block_linf_error, Real min_rhs_norm = 1.0e-6)
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

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_DIVERGENCE_FREE_NONZERO_RHS_MMS_HELPERS_H
