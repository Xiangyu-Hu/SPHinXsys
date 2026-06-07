#ifndef APHI_CONTACT_INTERFACE_DIAGNOSTIC_HELPERS_H
#define APHI_CONTACT_INTERFACE_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

/** Vol-weighted L2 mismatch ||lhs - rhs|| in a particle subset (uses names.true_residual = rhs - lhs). */
inline Real hostBlockMismatchNormInRegion(
    BaseParticles &particles, const AphiBlockNames &residual_block, const Vecd *positions, size_t total_real_particles,
    Real body_length, Real body_height, Real body_width, Real core_shell,
    const std::function<bool(const Vecd &)> &include_particle)
{
    return hostTrueResidualBlockNormInRegion(particles, residual_block, positions, total_real_particles, body_length,
                                             body_height, body_width, core_shell, include_particle);
}

struct AphiRegionalLhsRhsMismatch
{
    Real all_rel = 0.0;
    Real core_rel = 0.0;
    Real interface_band_rel = 0.0;
    Real core_away_interface_rel = 0.0;
    Real rhs_norm_all = 0.0;
    Real mismatch_norm_all = 0.0;
};

inline AphiRegionalLhsRhsMismatch measureBodyLhsRhsMismatch(
    SPHBody &body, Inner<> &inner, Contact<> *contact, const AphiVariableNames &names,
    const AphiLhsAssemblyOptions &options, Real body_length, Real body_height, Real body_width, Real core_shell,
    Real x_interface, Real interface_band_half_width, Real dp_0, const AphiBlockNames &input_block)
{
    AphiRegionalLhsRhsMismatch metrics;

    if (contact != nullptr)
    {
        AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply(
            body, inner, *contact, input_block, names.lhs, names.material, options.omega, options);
        apply.exec();
    }
    else
    {
        AphiApplyDynamicsBundle<MainExecutionPolicy> apply(body, inner, input_block, names.lhs, names.material,
                                                           options.omega, options);
        apply.exec();
    }

    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_mismatch(
        body, names.true_residual, names.rhs, names.lhs);
    compute_mismatch.exec();

    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    const auto all_particles = [](const Vecd &) { return true; };
    const auto core_particles = [&](const Vecd &position) {
        return isCoreParticle(position, body_length, body_height, body_width, core_shell);
    };
    const auto interface_band = [&](const Vecd &position) {
        return core_particles(position) && std::abs(position[0] - x_interface) <= interface_band_half_width;
    };
    const auto core_away_interface = [&](const Vecd &position) {
        return core_particles(position) && isAwayFromInterface(position, x_interface, dp_0);
    };

    metrics.mismatch_norm_all =
        hostBlockMismatchNormInRegion(particles, names.true_residual, positions, total_real_particles, body_length,
                                      body_height, body_width, core_shell, all_particles);
    metrics.rhs_norm_all = hostBlockNorm(particles, names.rhs, total_real_particles);
    metrics.all_rel = metrics.mismatch_norm_all / (metrics.rhs_norm_all + TinyReal);
    metrics.core_rel =
        hostBlockMismatchNormInRegion(particles, names.true_residual, positions, total_real_particles, body_length,
                                      body_height, body_width, core_shell, core_particles) /
        (metrics.rhs_norm_all + TinyReal);
    metrics.interface_band_rel =
        hostBlockMismatchNormInRegion(particles, names.true_residual, positions, total_real_particles, body_length,
                                      body_height, body_width, core_shell, interface_band) /
        (metrics.rhs_norm_all + TinyReal);
    metrics.core_away_interface_rel =
        hostBlockMismatchNormInRegion(particles, names.true_residual, positions, total_real_particles, body_length,
                                      body_height, body_width, core_shell, core_away_interface) /
        (metrics.rhs_norm_all + TinyReal);
    return metrics;
}

inline void assembleTwoBodyCoupledRhsFromExact(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                               const AphiLhsAssemblyOptions &options)
{
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left_exact(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right_exact(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_exact_reference(case_setup.left_body, names.r_hat,
                                                                                names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_exact_reference(case_setup.right_body, names.r_hat,
                                                                                 names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_lhs_to_rhs(case_setup.left_body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_lhs_to_rhs(case_setup.right_body, names.rhs,
                                                                            names.lhs);

    apply_left_exact.exec();
    apply_right_exact.exec();
    copy_left_exact_reference.exec();
    copy_right_exact_reference.exec();
    copy_left_lhs_to_rhs.exec();
    copy_right_lhs_to_rhs.exec();
    case_setup.updateRelations();
}

inline void zeroLeftSolutionOnly(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names)
{
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left(case_setup.left_body, names.solution);
    zero_left.exec();
    syncAphiBlockToDevice(case_setup.left_body.getBaseParticles(), names.solution);
}

inline void pinRightSolutionFromRhat(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names)
{
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> pin_right(case_setup.right_body, names.r_hat, names.solution);
    pin_right.exec();
    syncAphiBlockToDevice(case_setup.right_body.getBaseParticles(), names.solution);
}

inline void syncTwoBodySolutionToDevice(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names)
{
    syncAphiBlockToDevice(case_setup.left_body.getBaseParticles(), names.solution);
    syncAphiBlockToDevice(case_setup.right_body.getBaseParticles(), names.solution);
}

inline void pinTwoBodySolutionFromRhat(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names)
{
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> pin_left(case_setup.left_body, names.r_hat, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> pin_right(case_setup.right_body, names.r_hat, names.solution);
    pin_left.exec();
    pin_right.exec();
    syncTwoBodySolutionToDevice(case_setup, names);
}

inline void zeroTwoBodySolution(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names)
{
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left(case_setup.left_body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right(case_setup.right_body, names.solution);
    zero_left.exec();
    zero_right.exec();
    syncTwoBodySolutionToDevice(case_setup, names);
}

inline void assembleMonolithicRhsFromExact(AphiLhsTestBody &test_body, const AphiVariableNames &names,
                                           const AphiLhsAssemblyOptions &options)
{
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_exact(test_body.body, test_body.inner(), names.solution,
                                                             names.lhs, names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact_reference(test_body.body, names.r_hat, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);

    apply_exact.exec();
    copy_exact_reference.exec();
    copy_lhs_to_rhs.exec();
    test_body.updateRelations();
}

struct AphiSplitMonoLhsComparison
{
    size_t matched_particles = 0;
    size_t missing_particles = 0;
    Real max_abs_diff = 0.0;
    Real interface_band_max_abs_diff = 0.0;
    size_t interface_band_matched = 0;
};

inline void pinMonolithicSolutionFromRhat(AphiLhsTestBody &test_body, const AphiVariableNames &names)
{
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> pin_solution(test_body.body, names.r_hat, names.solution);
    pin_solution.exec();
}

inline Real hostMaxAbsBlockDifference(BaseParticles &particles, const AphiBlockNames &lhs_block,
                                      const AphiBlockNames &rhs_block, size_t total_real_particles)
{
    syncAphiBlockToHost(particles, lhs_block);
    syncAphiBlockToHost(particles, rhs_block);
    const Vecd *lhs_a_real = particles.getVariableDataByName<Vecd>(lhs_block.a_real);
    const Vecd *lhs_a_imag = particles.getVariableDataByName<Vecd>(lhs_block.a_imag);
    const Real *lhs_phi_real = particles.getVariableDataByName<Real>(lhs_block.phi_real);
    const Real *lhs_phi_imag = particles.getVariableDataByName<Real>(lhs_block.phi_imag);
    const Vecd *rhs_a_real = particles.getVariableDataByName<Vecd>(rhs_block.a_real);
    const Vecd *rhs_a_imag = particles.getVariableDataByName<Vecd>(rhs_block.a_imag);
    const Real *rhs_phi_real = particles.getVariableDataByName<Real>(rhs_block.phi_real);
    const Real *rhs_phi_imag = particles.getVariableDataByName<Real>(rhs_block.phi_imag);
    Real max_diff = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        max_diff = std::max(max_diff, (lhs_a_real[i] - rhs_a_real[i]).norm());
        max_diff = std::max(max_diff, (lhs_a_imag[i] - rhs_a_imag[i]).norm());
        max_diff = std::max(max_diff, std::abs(lhs_phi_real[i] - rhs_phi_real[i]));
        max_diff = std::max(max_diff, std::abs(lhs_phi_imag[i] - rhs_phi_imag[i]));
    }
    return max_diff;
}

inline Real twoBodyMaxAbsSolutionRhatDifference(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names)
{
    const size_t n_left = case_setup.left_body.getBaseParticles().TotalRealParticles();
    const size_t n_right = case_setup.right_body.getBaseParticles().TotalRealParticles();
    return std::max(hostMaxAbsBlockDifference(case_setup.left_body.getBaseParticles(), names.solution, names.r_hat,
                                              n_left),
                    hostMaxAbsBlockDifference(case_setup.right_body.getBaseParticles(), names.solution, names.r_hat,
                                              n_right));
}

inline AphiSplitMonoLhsComparison compareSplitContactLhsToMonolithic(
    AphiLhsTestBody &mono_body, AphiTwoBodyInterfaceCase &split_case, const AphiVariableNames &names,
    const AphiLhsAssemblyOptions &options, Real body_length, Real body_height, Real body_width, Real core_shell,
    Real x_interface, Real interface_band_half_width, Real dp_0)
{
    AphiSplitMonoLhsComparison summary;

    pinMonolithicSolutionFromRhat(mono_body, names);
    pinTwoBodySolutionFromRhat(split_case, names);
    mono_body.updateRelations();
    split_case.updateRelations();

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_mono(mono_body.body, mono_body.inner(), names.solution, names.lhs,
                                                            names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left(
        split_case.left_body, split_case.left_inner(), split_case.left_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right(
        split_case.right_body, split_case.right_inner(), split_case.right_contact(), names.solution, names.lhs,
        names.material, options.omega, options);

    apply_mono.exec();
    apply_left.exec();
    apply_right.exec();
    mono_body.updateRelations();
    split_case.updateRelations();

    AphiBlockMapByPosition mono_lhs;
    AphiBlockMapByPosition split_lhs;
    collectCoreBlockByPosition(mono_lhs, mono_body.body.getBaseParticles(), names.lhs, body_length, body_height,
                               body_width, core_shell, x_interface, dp_0);
    collectCoreBlockByPosition(split_lhs, split_case.left_body.getBaseParticles(), names.lhs, body_length, body_height,
                               body_width, core_shell, x_interface, dp_0);
    collectCoreBlockByPosition(split_lhs, split_case.right_body.getBaseParticles(), names.lhs, body_length, body_height,
                               body_width, core_shell, x_interface, dp_0);

    summary.max_abs_diff =
        maxAbsBlockDifference(mono_lhs, split_lhs, summary.matched_particles, summary.missing_particles);

    for (const auto &entry : mono_lhs)
    {
        const auto it = split_lhs.find(entry.first);
        if (it == split_lhs.end())
        {
            continue;
        }
        if (std::abs(entry.first.position[0] - x_interface) > interface_band_half_width)
        {
            continue;
        }
        const AphiContactBlockByPosition &mono_block = entry.second;
        const AphiContactBlockByPosition &split_block = it->second;
        summary.interface_band_matched += 1;
        summary.interface_band_max_abs_diff =
            std::max(summary.interface_band_max_abs_diff, (mono_block.a_real - split_block.a_real).norm());
        summary.interface_band_max_abs_diff =
            std::max(summary.interface_band_max_abs_diff, (mono_block.a_imag - split_block.a_imag).norm());
        summary.interface_band_max_abs_diff =
            std::max(summary.interface_band_max_abs_diff, std::abs(mono_block.phi_real - split_block.phi_real));
        summary.interface_band_max_abs_diff =
            std::max(summary.interface_band_max_abs_diff, std::abs(mono_block.phi_imag - split_block.phi_imag));
    }

    return summary;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_INTERFACE_DIAGNOSTIC_HELPERS_H
