#ifndef APHI_CONTACT_WORKSPACE_CONTAMINATION_HELPERS_H
#define APHI_CONTACT_WORKSPACE_CONTAMINATION_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline Real hostRelativeBlockDifference(SPHBody &body, const AphiBlockNames &reference_block,
                                        const AphiBlockNames &other_block)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    syncAphiBlockToHost(particles, reference_block);
    syncAphiBlockToHost(particles, other_block);
    Real sum_squared_diff = 0.0;
    Real sum_squared_ref = 0.0;
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *ref_a_real = particles.getVariableDataByName<Vecd>(reference_block.a_real);
    const Vecd *ref_a_imag = particles.getVariableDataByName<Vecd>(reference_block.a_imag);
    const Real *ref_phi_real = particles.getVariableDataByName<Real>(reference_block.phi_real);
    const Real *ref_phi_imag = particles.getVariableDataByName<Real>(reference_block.phi_imag);
    const Vecd *other_a_real = particles.getVariableDataByName<Vecd>(other_block.a_real);
    const Vecd *other_a_imag = particles.getVariableDataByName<Vecd>(other_block.a_imag);
    const Real *other_phi_real = particles.getVariableDataByName<Real>(other_block.phi_real);
    const Real *other_phi_imag = particles.getVariableDataByName<Real>(other_block.phi_imag);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        const Vecd da_real = ref_a_real[i] - other_a_real[i];
        const Vecd da_imag = ref_a_imag[i] - other_a_imag[i];
        const Real dphi_real = ref_phi_real[i] - other_phi_real[i];
        const Real dphi_imag = ref_phi_imag[i] - other_phi_imag[i];
        sum_squared_diff += vol[i] * (da_real.squaredNorm() + da_imag.squaredNorm() + dphi_real * dphi_real +
                                      dphi_imag * dphi_imag);
        sum_squared_ref += vol[i] * (ref_a_real[i].squaredNorm() + ref_a_imag[i].squaredNorm() +
                                     ref_phi_real[i] * ref_phi_real[i] + ref_phi_imag[i] * ref_phi_imag[i]);
    }
    return std::sqrt(sum_squared_diff) / (std::sqrt(sum_squared_ref) + TinyReal);
}

inline void applyZToSearchOnLeftBody(AphiTwoBodyInterfaceCase &case_setup, const AphiBlockNames &z_block,
                                     const AphiBlockNames &output_block, const AphiMaterialNames &material_names,
                                     const AphiLhsAssemblyOptions &options, bool use_block_diagonal)
{
    if (use_block_diagonal)
    {
        AphiApplyContactBlockDiagonalDynamicsBundle<MainExecutionPolicy> apply(
            case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), z_block, output_block,
            material_names, options.omega, options);
        apply.exec();
    }
    else
    {
        AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply(
            case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), z_block, output_block,
            material_names, options.omega, options);
        apply.exec();
    }
}

struct WorkspaceContaminationMetrics
{
    Real full_apply_rel_diff = 0.0;
    Real block_diagonal_rel_diff = 0.0;
    Real reference_output_norm = 0.0;
};

inline WorkspaceContaminationMetrics runWorkspaceContaminationCheck(AphiTwoBodyInterfaceCase &case_setup,
                                                                    const AphiVariableNames &names,
                                                                    const AphiLhsAssemblyOptions &options)
{
    WorkspaceContaminationMetrics metrics;
    const AphiBlockNames z0 = aphiGMRESZBlock(0);
    AphiVariableNames z0_field_names = names;
    z0_field_names.solution = z0;

    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_left_z0(case_setup.left_body, z0_field_names);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right_z0(case_setup.right_body, z0);
    assign_left_z0.exec();
    zero_right_z0.exec();
    case_setup.updateRelations();

    applyZToSearchOnLeftBody(case_setup, z0, names.search, names.material, options, false);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> save_full_reference(case_setup.left_body, names.t, names.search);
    save_full_reference.exec();
    metrics.reference_output_norm =
        hostBlockNorm(case_setup.left_body.getBaseParticles(), names.t,
                      case_setup.left_body.getBaseParticles().TotalRealParticles());

    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_dirty_right_z0(case_setup.right_body,
                                                                                          z0_field_names);
    assign_dirty_right_z0.exec();
    case_setup.updateRelations();

    applyZToSearchOnLeftBody(case_setup, z0, names.search, names.material, options, false);
    metrics.full_apply_rel_diff = hostRelativeBlockDifference(case_setup.left_body, names.t, names.search);

    applyZToSearchOnLeftBody(case_setup, z0, names.search, names.material, options, true);
    metrics.block_diagonal_rel_diff = hostRelativeBlockDifference(case_setup.left_body, names.t, names.search);
    return metrics;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_WORKSPACE_CONTAMINATION_HELPERS_H
