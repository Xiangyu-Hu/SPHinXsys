#ifndef APHI_CONTACT_INTERFACE_PENALTY_BREAKDOWN_HELPERS_H
#define APHI_CONTACT_INTERFACE_PENALTY_BREAKDOWN_HELPERS_H

#include "electromagnetic_dynamics/diagnostics/aphi_contact_interface_mms_consistency_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactInterfacePenaltyBreakdownRow
{
    Real eta_a = 0.0;
    Real lambda_a = 0.0;
    AphiContactRegionalMmsDefectMetrics interface_base{};
    AphiContactRegionalMmsDefectMetrics interface_full{};
    AphiContactRegionalMmsDefectMetrics interface_penalty_only{};
    AphiContactRegionalMmsDefectMetrics boundary_penalty_only{};
};

inline void execTwoBodyContactApply(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                    const AphiLhsAssemblyOptions &options, Real omega)
{
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, omega, options);
    apply_left.exec();
    apply_right.exec();
    case_setup.updateRelations();
}

inline AphiContactRegionalMmsDefectMetrics hostTwoBodyRegionalLhsRhsDefectWithBlock(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names, const AphiBlockNames &lhs_block,
    Real body_length, Real body_height, Real body_width, Real core_shell,
    const std::function<bool(const Vecd &, bool is_core)> &in_region)
{
    AphiContactRegionalMmsDefectMetrics metrics;
    Real diff_squared = 0.0;
    Real ref_squared = 0.0;

    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        syncAphiBlockToHost(particles, lhs_block);
        syncAphiBlockToHost(particles, names.rhs);
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const Vecd *lhs_a_real = particles.getVariableDataByName<Vecd>(lhs_block.a_real);
        const Vecd *lhs_a_imag = particles.getVariableDataByName<Vecd>(lhs_block.a_imag);
        const Real *lhs_phi_real = particles.getVariableDataByName<Real>(lhs_block.phi_real);
        const Real *lhs_phi_imag = particles.getVariableDataByName<Real>(lhs_block.phi_imag);
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

inline AphiContactRegionalMmsDefectMetrics hostTwoBodyRegionalLhsRhsDefect(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names, Real body_length, Real body_height,
    Real body_width, Real core_shell, const std::function<bool(const Vecd &, bool is_core)> &in_region)
{
    return hostTwoBodyRegionalLhsRhsDefectWithBlock(case_setup, names, names.lhs, body_length, body_height, body_width,
                                                    core_shell, in_region);
}

inline AphiContactRegionalMmsDefectMetrics hostTwoBodyRegionalLhsDifference(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names, const AphiBlockNames &lhs_a,
    const AphiBlockNames &lhs_b, Real body_length, Real body_height, Real body_width, Real core_shell,
    const std::function<bool(const Vecd &, bool is_core)> &in_region)
{
    AphiContactRegionalMmsDefectMetrics metrics;
    Real diff_squared = 0.0;
    Real ref_squared = 0.0;

    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        syncAphiBlockToHost(particles, lhs_a);
        syncAphiBlockToHost(particles, lhs_b);
        syncAphiBlockToHost(particles, names.rhs);
        syncVariableToHost<Vecd>(particles, "Position");
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
        const Vecd *a_a_real = particles.getVariableDataByName<Vecd>(lhs_a.a_real);
        const Vecd *a_a_imag = particles.getVariableDataByName<Vecd>(lhs_a.a_imag);
        const Real *a_phi_real = particles.getVariableDataByName<Real>(lhs_a.phi_real);
        const Real *a_phi_imag = particles.getVariableDataByName<Real>(lhs_a.phi_imag);
        const Vecd *b_a_real = particles.getVariableDataByName<Vecd>(lhs_b.a_real);
        const Vecd *b_a_imag = particles.getVariableDataByName<Vecd>(lhs_b.a_imag);
        const Real *b_phi_real = particles.getVariableDataByName<Real>(lhs_b.phi_real);
        const Real *b_phi_imag = particles.getVariableDataByName<Real>(lhs_b.phi_imag);
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
            const Vecd da_re = a_a_real[i] - b_a_real[i];
            const Vecd da_im = a_a_imag[i] - b_a_imag[i];
            const Real dphi_re = a_phi_real[i] - b_phi_real[i];
            const Real dphi_im = a_phi_imag[i] - b_phi_imag[i];
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

inline AphiContactInterfacePenaltyBreakdownRow runContactInterfacePenaltyBreakdownRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma, Real nu, Real omega, Real phi_gauge_penalty, Real eta_a,
    const AphiCoreOperatorScaleMetrics &scale_metrics, Real x_interface, Real interface_band_half_width)
{
    AphiContactInterfacePenaltyBreakdownRow row;
    row.eta_a = eta_a;
    row.lambda_a = eta_a > TinyReal ? lambdaAFromEtaA(eta_a, scale_metrics) : 0.0;

    AphiTwoBodyInterfaceCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    setupTwoBodyDivFreeMmsFields(case_setup, names, sigma, nu, AphiDivFreeValidationFieldKind::Az2D);
    case_setup.updateRelations();

    AphiLhsAssemblyOptions base_options;
    base_options.omega = omega;
    base_options.use_phi_gauge_penalty = true;
    base_options.phi_gauge_penalty = phi_gauge_penalty;
    base_options.use_a_divergence_penalty = false;

    AphiLhsAssemblyOptions full_options = base_options;
    full_options.use_a_divergence_penalty = eta_a > TinyReal;
    full_options.a_divergence_penalty = row.lambda_a;
    full_options.a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;

    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_lhs_to_rhs(case_setup.left_body, names.rhs,
                                                                               names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_lhs_to_rhs(case_setup.right_body, names.rhs,
                                                                              names.lhs);

    execTwoBodyContactApply(case_setup, names, base_options, omega);
    copy_left_lhs_to_rhs.exec();
    copy_right_lhs_to_rhs.exec();
    case_setup.updateRelations();
    execTwoBodyContactApply(case_setup, names, base_options, omega);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_base_lhs(case_setup.left_body, names.t, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_base_lhs(case_setup.right_body, names.t, names.lhs);
    copy_left_base_lhs.exec();
    copy_right_base_lhs.exec();

    execTwoBodyContactApply(case_setup, names, full_options, omega);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_full_lhs(case_setup.left_body, names.v, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_full_lhs(case_setup.right_body, names.v, names.lhs);
    copy_left_full_lhs.exec();
    copy_right_full_lhs.exec();
    case_setup.updateRelations();

    const auto interface_band = [&](const Vecd &position, bool is_core) {
        return is_core && std::abs(position[0] - x_interface) <= interface_band_half_width + TinyReal;
    };
    const auto boundary = [&](const Vecd &position, bool is_core) { return !is_core; };

    row.interface_base = hostTwoBodyRegionalLhsRhsDefectWithBlock(
        case_setup, names, names.t, body_length, body_height, body_width, core_shell, interface_band);
    row.interface_full = hostTwoBodyRegionalLhsRhsDefectWithBlock(
        case_setup, names, names.v, body_length, body_height, body_width, core_shell, interface_band);
    row.interface_penalty_only = hostTwoBodyRegionalLhsDifference(
        case_setup, names, names.v, names.t, body_length, body_height, body_width, core_shell, interface_band);
    row.boundary_penalty_only = hostTwoBodyRegionalLhsDifference(
        case_setup, names, names.v, names.t, body_length, body_height, body_width, core_shell, boundary);
    return row;
}

inline void printContactInterfacePenaltyBreakdownRow(const char *test_name,
                                                     const AphiContactInterfacePenaltyBreakdownRow &row)
{
    std::cout << test_name << " eta_A=" << row.eta_a << " lambda_A=" << row.lambda_a
              << " interface_base_rel=" << row.interface_base.relative_l2
              << " interface_full_rel=" << row.interface_full.relative_l2
              << " interface_penalty_only_rel=" << row.interface_penalty_only.relative_l2
              << " interface_penalty_only_max=" << row.interface_penalty_only.max_abs_block_diff
              << " boundary_penalty_only_rel=" << row.boundary_penalty_only.relative_l2
              << " boundary_penalty_only_max=" << row.boundary_penalty_only.max_abs_block_diff << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_INTERFACE_PENALTY_BREAKDOWN_HELPERS_H
