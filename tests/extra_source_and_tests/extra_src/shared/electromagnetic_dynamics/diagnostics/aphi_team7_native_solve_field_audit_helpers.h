#ifndef APHI_TEAM7_NATIVE_SOLVE_FIELD_AUDIT_HELPERS_H
#define APHI_TEAM7_NATIVE_SOLVE_FIELD_AUDIT_HELPERS_H

#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_team7_native_reload_geometry_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiTeam7NativeAirSolveFieldAudit
{
    bool air_multiresolution = false;
    size_t air_particles = 0;
    Real reference_smoothing_length_m = 0.0;
    Real h_ratio_min = 1.0;
    Real h_ratio_max = 1.0;
    Real h_ratio_mean = 1.0;
    Real solution_block_l2 = 0.0;
    Real rhs_block_l2 = 0.0;
    Real lhs_block_l2 = 0.0;
    Real air_only_residual_l2 = 0.0;
    Real air_only_relative_residual = 0.0;
    Real max_abs_a_real = 0.0;
    Real max_abs_a_imag = 0.0;
    Real max_abs_phi_real = 0.0;
    UnsignedInt neighbor_count_min = 0;
    UnsignedInt neighbor_count_max = 0;
    Real neighbor_count_mean = 0.0;
};

inline void team7NativeApplyAirLhsOnCase(AphiTeam7NativeReloadContactCase &case_setup, const AphiVariableNames &names,
                                         const AphiLhsAssemblyOptions &options)
{
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_lhs(case_setup.air_body(), names.lhs);
    InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<Team7NativeAdaptiveAirInner>> apply_inner(DynamicsArgs(
        case_setup.air_inner(), names.solution, names.lhs, names.material, options.omega, options, Real(0.01)));
    InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<Team7AirToCoilContact>> apply_coil(DynamicsArgs(
        case_setup.air_to_coil_contact(), names.solution, names.lhs, names.material, options.omega, options, Real(0.01)));
    InteractionDynamicsCK<MainExecutionPolicy, AphiApplyCK<Team7AirToPlateContact>> apply_plate(DynamicsArgs(
        case_setup.air_to_plate_contact(), names.solution, names.lhs, names.material, options.omega, options, Real(0.01)));
    zero_lhs.exec();
    apply_inner.exec();
    apply_coil.exec();
    apply_plate.exec();
}

inline AphiTeam7NativeAirSolveFieldAudit auditTeam7NativeAirSolveFields(AphiTeam7NativeReloadContactCase &case_setup,
                                                                        const AphiVariableNames &names,
                                                                        const AphiLhsAssemblyOptions &options)
{
    AphiTeam7NativeAirSolveFieldAudit audit;
    audit.air_multiresolution = case_setup.geometry_config.local_refinement_levels > 0;
    BaseParticles &particles = case_setup.air_body().getBaseParticles();
    audit.air_particles = particles.TotalRealParticles();
    audit.reference_smoothing_length_m = case_setup.air_body().getSPHAdaptation().ReferenceSmoothingLength();

    syncVariableToHost<Vecd>(particles, "Position");
    const size_t count = particles.TotalRealParticles();
    const auto all_particles = [](const Vecd &) { return true; };
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");

    audit.solution_block_l2 =
        hostParticleRegionBlockNorm(particles, names.solution, positions, count, all_particles);
    audit.rhs_block_l2 = hostParticleRegionBlockNorm(particles, names.rhs, positions, count, all_particles);

    team7NativeApplyAirLhsOnCase(case_setup, names, options);
    audit.lhs_block_l2 = hostParticleRegionBlockNorm(particles, names.lhs, positions, count, all_particles);

    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> air_residual(
        case_setup.air_body(), names.true_residual, names.rhs, names.lhs);
    air_residual.exec();
    audit.air_only_residual_l2 =
        hostParticleRegionBlockNorm(particles, names.true_residual, positions, count, all_particles);
    audit.air_only_relative_residual = audit.air_only_residual_l2 / (audit.rhs_block_l2 + TinyReal);

    syncAphiBlockToHost(particles, names.solution);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(names.solution.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(names.solution.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(names.solution.phi_real);
    for (size_t i = 0; i != count; ++i)
    {
        audit.max_abs_a_real = std::max(audit.max_abs_a_real, a_real[i].norm());
        audit.max_abs_a_imag = std::max(audit.max_abs_a_imag, a_imag[i].norm());
        audit.max_abs_phi_real = std::max(audit.max_abs_phi_real, std::abs(phi_real[i]));
    }

    syncVariableToHost<UnsignedInt>(particles, "NeighborSize");
    const UnsignedInt *neighbor_size = particles.getVariableDataByName<UnsignedInt>("NeighborSize");
    UnsignedInt neighbor_sum = 0;
    audit.neighbor_count_min = std::numeric_limits<UnsignedInt>::max();
    audit.neighbor_count_max = 0;
    for (size_t i = 0; i != count; ++i)
    {
        const UnsignedInt n_neighbors = neighbor_size[i];
        neighbor_sum += n_neighbors;
        audit.neighbor_count_min = std::min(audit.neighbor_count_min, n_neighbors);
        audit.neighbor_count_max = std::max(audit.neighbor_count_max, n_neighbors);
    }
    audit.neighbor_count_mean = static_cast<Real>(neighbor_sum) / static_cast<Real>(count);

    syncVariableToHost<Real>(particles, "SmoothingLengthRatio");
    const Real *h_ratio = particles.getVariableDataByName<Real>("SmoothingLengthRatio");
    audit.h_ratio_min = h_ratio[0];
    audit.h_ratio_max = h_ratio[0];
    Real h_sum = 0.0;
    for (size_t i = 0; i != count; ++i)
    {
        audit.h_ratio_min = std::min(audit.h_ratio_min, h_ratio[i]);
        audit.h_ratio_max = std::max(audit.h_ratio_max, h_ratio[i]);
        h_sum += h_ratio[i];
    }
    audit.h_ratio_mean = h_sum / static_cast<Real>(count);
    return audit;
}

inline void printTeam7NativeAirSolveFieldAudit(const char *label, const AphiTeam7NativeAirSolveFieldAudit &audit)
{
    std::cout << "team7_native_air_solve_audit label=" << label << " air_mr=" << (audit.air_multiresolution ? 1 : 0)
              << " air_particles=" << audit.air_particles
              << " h_ref_m=" << audit.reference_smoothing_length_m << " h_ratio_min=" << audit.h_ratio_min
              << " h_ratio_max=" << audit.h_ratio_max << " h_ratio_mean=" << audit.h_ratio_mean
              << " neighbor_min=" << audit.neighbor_count_min << " neighbor_max=" << audit.neighbor_count_max
              << " neighbor_mean=" << audit.neighbor_count_mean << " solution_l2=" << audit.solution_block_l2
              << " rhs_l2=" << audit.rhs_block_l2 << " lhs_l2=" << audit.lhs_block_l2
              << " air_residual_l2=" << audit.air_only_residual_l2
              << " air_rel_res=" << audit.air_only_relative_residual << " max_abs_A_real=" << audit.max_abs_a_real
              << " max_abs_A_imag=" << audit.max_abs_a_imag << " max_abs_phi_real=" << audit.max_abs_phi_real
              << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_TEAM7_NATIVE_SOLVE_FIELD_AUDIT_HELPERS_H
