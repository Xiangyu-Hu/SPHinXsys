#ifndef APHI_CONTACT_INTERFACE_MMS_CONSISTENCY_HELPERS_H
#define APHI_CONTACT_INTERFACE_MMS_CONSISTENCY_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_contact_a_divergence_penalty_two_body_mms_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactRegionalMmsDefectMetrics
{
    size_t particle_count = 0;
    Real max_abs_block_diff = 0.0;
    Real weighted_l2_diff = 0.0;
    Real weighted_ref_l2 = 0.0;
    Real relative_l2 = 0.0;
};

struct AphiContactInterfaceMmsConsistencyRow
{
    Real eta_a = 0.0;
    Real lambda_a = 0.0;
    Real global_relative_l2 = 0.0;
    AphiContactRegionalMmsDefectMetrics stencil_safe_core{};
    AphiContactRegionalMmsDefectMetrics interface_band{};
    AphiContactRegionalMmsDefectMetrics boundary{};
};

inline AphiContactRegionalMmsDefectMetrics hostTwoBodyRegionalMmsDefectMetrics(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names, Real body_length, Real body_height,
    Real body_width, Real core_shell, Real x_interface, Real interface_band_half_width,
    const std::function<bool(const Vecd &, bool is_core)> &in_region)
{
    AphiContactRegionalMmsDefectMetrics metrics;
    Real diff_squared = 0.0;
    Real ref_squared = 0.0;

    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        syncAphiBlockToHost(particles, names.lhs);
        syncAphiBlockToHost(particles, names.rhs);
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const Vecd *lhs_a_real = particles.getVariableDataByName<Vecd>(names.lhs.a_real);
        const Vecd *lhs_a_imag = particles.getVariableDataByName<Vecd>(names.lhs.a_imag);
        const Real *lhs_phi_real = particles.getVariableDataByName<Real>(names.lhs.phi_real);
        const Real *lhs_phi_imag = particles.getVariableDataByName<Real>(names.lhs.phi_imag);
        const Vecd *rhs_a_real = particles.getVariableDataByName<Vecd>(names.rhs.a_real);
        const Vecd *rhs_a_imag = particles.getVariableDataByName<Vecd>(names.rhs.a_imag);
        const Real *rhs_phi_real = particles.getVariableDataByName<Real>(names.rhs.phi_real);
        const Real *rhs_phi_imag = particles.getVariableDataByName<Real>(names.rhs.phi_imag);

        for (size_t i = 0; i != total_real_particles; ++i)
        {
            const bool is_core =
                isCoreParticle(positions[i], body_length, body_height, body_width, core_shell);
            if (!in_region(positions[i], is_core))
            {
                continue;
            }
            const Vecd da_re = lhs_a_real[i] - rhs_a_real[i];
            const Vecd da_im = lhs_a_imag[i] - rhs_a_imag[i];
            const Real dphi_re = lhs_phi_real[i] - rhs_phi_real[i];
            const Real dphi_im = lhs_phi_imag[i] - rhs_phi_imag[i];
            const Real block_diff = std::sqrt(da_re.squaredNorm() + da_im.squaredNorm() + dphi_re * dphi_re +
                                              dphi_im * dphi_im);
            const Real block_ref = std::sqrt(rhs_a_real[i].squaredNorm() + rhs_a_imag[i].squaredNorm() +
                                             rhs_phi_real[i] * rhs_phi_real[i] + rhs_phi_imag[i] * rhs_phi_imag[i]);
            metrics.particle_count += 1;
            metrics.max_abs_block_diff = std::max(metrics.max_abs_block_diff, block_diff);
            diff_squared += vol[i] * block_diff * block_diff;
            ref_squared += vol[i] * block_ref * block_ref;
        }
    }

    metrics.weighted_l2_diff = std::sqrt(diff_squared);
    metrics.weighted_ref_l2 = std::sqrt(ref_squared);
    metrics.relative_l2 = metrics.weighted_l2_diff / (metrics.weighted_ref_l2 + TinyReal);
    return metrics;
}

inline AphiContactInterfaceMmsConsistencyRow runContactInterfaceMmsConsistencyRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real x_interface, Real interface_band_half_width,
    AphiDivFreeValidationFieldKind field_kind,
    AphiContactADivergencePenaltyStencilMode penalty_stencil =
        AphiContactADivergencePenaltyStencilMode::InnerOnly)
{
    AphiContactInterfaceMmsConsistencyRow row;
    row.eta_a = eta_a;
    row.lambda_a = eta_a > TinyReal ? lambdaAFromEtaA(eta_a, scale_metrics) : 0.0;

    AphiTwoBodyInterfaceCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    setupTwoBodyDivFreeMmsFields(case_setup, names, sigma, nu, field_kind);
    case_setup.updateRelations();

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = eta_a > TinyReal;
    options.a_divergence_penalty = row.lambda_a;
    options.a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;
    options.contact_a_divergence_penalty_stencil = penalty_stencil;

    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, omega, options);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_lhs_to_rhs(case_setup.left_body, names.rhs,
                                                                               names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_lhs_to_rhs(case_setup.right_body, names.rhs,
                                                                              names.lhs);

    apply_left.exec();
    apply_right.exec();
    copy_left_lhs_to_rhs.exec();
    copy_right_lhs_to_rhs.exec();
    case_setup.updateRelations();
    apply_left.exec();
    apply_right.exec();
    case_setup.updateRelations();

    const auto stencil_safe_core = [&](const Vecd &position, bool is_core) {
        return is_core && isGradStencilSafeFromInterface(position, x_interface, dp_0);
    };
    const auto interface_band = [&](const Vecd &position, bool is_core) {
        return is_core && std::abs(position[0] - x_interface) <= interface_band_half_width + TinyReal;
    };
    const auto boundary = [&](const Vecd &position, bool is_core) { return !is_core; };

    row.stencil_safe_core = hostTwoBodyRegionalMmsDefectMetrics(
        case_setup, names, body_length, body_height, body_width, core_shell, x_interface, interface_band_half_width,
        stencil_safe_core);
    row.interface_band = hostTwoBodyRegionalMmsDefectMetrics(
        case_setup, names, body_length, body_height, body_width, core_shell, x_interface, interface_band_half_width,
        interface_band);
    row.boundary = hostTwoBodyRegionalMmsDefectMetrics(
        case_setup, names, body_length, body_height, body_width, core_shell, x_interface, interface_band_half_width,
        boundary);

    Real global_diff_squared = row.stencil_safe_core.weighted_l2_diff * row.stencil_safe_core.weighted_l2_diff +
                               row.interface_band.weighted_l2_diff * row.interface_band.weighted_l2_diff +
                               row.boundary.weighted_l2_diff * row.boundary.weighted_l2_diff;
    Real global_ref_squared = row.stencil_safe_core.weighted_ref_l2 * row.stencil_safe_core.weighted_ref_l2 +
                              row.interface_band.weighted_ref_l2 * row.interface_band.weighted_ref_l2 +
                              row.boundary.weighted_ref_l2 * row.boundary.weighted_ref_l2;
    row.global_relative_l2 = std::sqrt(global_diff_squared) / (std::sqrt(global_ref_squared) + TinyReal);
    return row;
}

inline void printContactInterfaceMmsConsistencyRow(const char *test_name, const AphiContactInterfaceMmsConsistencyRow &row)
{
    std::cout << test_name << " eta_A=" << row.eta_a << " lambda_A=" << row.lambda_a
              << " global_rel_l2=" << row.global_relative_l2 << " core_safe_rel_l2=" << row.stencil_safe_core.relative_l2
              << " core_safe_max=" << row.stencil_safe_core.max_abs_block_diff
              << " core_safe_N=" << row.stencil_safe_core.particle_count
              << " interface_rel_l2=" << row.interface_band.relative_l2
              << " interface_max=" << row.interface_band.max_abs_block_diff
              << " interface_N=" << row.interface_band.particle_count
              << " boundary_rel_l2=" << row.boundary.relative_l2 << " boundary_max=" << row.boundary.max_abs_block_diff
              << " boundary_N=" << row.boundary.particle_count << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_INTERFACE_MMS_CONSISTENCY_HELPERS_H
