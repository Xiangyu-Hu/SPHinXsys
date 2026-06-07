/**
 * Stage 10.5 diagnostic: left conductor field error after coupled Contact GMRES.
 * 10.5-A gauge/mean-offset; 10.5-C core vs interface band; 10.5-E restricted monolithic.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_contact_left_field_error_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

inline void runMonolithicInterfaceMms(AphiLhsTestBody &test_body, const AphiVariableNames &names,
                                      const AphiLhsAssemblyOptions &options,
                                      const AphiMatrixFreeSolverOptions &solver_options)
{
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_exact_reference(test_body.body, names.r_hat, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_exact(test_body.body, test_body.inner(), names.solution,
                                                             names.lhs, names.material, options.omega, options);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                      solver_options);

    apply_exact.exec();
    copy_exact_reference.exec();
    copy_lhs_to_rhs.exec();
    zero_solution.exec();
    test_body.updateRelations();
    (void)solver.solve();
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real sigma_conductor = 2.0;
    const Real sigma_air = 1.0e-4;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real tolerance = 1.0e-5;
    const Real x_interface = 0.5;
    const Real interface_band_half_width = 2.0 * dp_0;

    AphiVariableNames names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    const AphiMatrixFreeSolverOptions mono_solver_options = defaultMatrixFreeGMRESOptions(tolerance, 50, 30);
    const AphiMatrixFreeSolverOptions contact_solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 100);

    AphiLhsTestBody mono_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_mono(mono_body.body, sigma_conductor, nu,
                                                                                  names);
    StateDynamics<MainExecutionPolicy, AssignPiecewiseSigmaHalfSpaceCK> assign_mono_material(
        mono_body.body, x_interface, sigma_conductor, sigma_air, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignInterfaceFluxMatchedPhiFieldsCK> assign_mono_exact(
        mono_body.body, x_interface, sigma_conductor, sigma_air, names.solution);
    initialize_mono.exec();
    assign_mono_material.exec();
    assign_mono_exact.exec();
    mono_body.updateRelations();
    runMonolithicInterfaceMms(mono_body, names, options, mono_solver_options);

    AphiTwoBodyInterfaceCase contact_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    setupTwoBodyInterfaceMmsFields(contact_case, names, sigma_conductor, sigma_air, nu, x_interface);
    contact_case.updateRelations();

    const AphiInterfaceMmsRunMetrics contact_metrics = runTwoBodyCoupledContactInterfaceMms(
        contact_case, names, options, contact_solver_options, body_length, body_height, body_width, core_shell,
        x_interface, interface_band_half_width);

    BaseParticles &left_particles = contact_case.left_body.getBaseParticles();
    syncVariableToHost<Vecd>(left_particles, "Position");
    const size_t left_total_real_particles = left_particles.TotalRealParticles();
    const Vecd *left_positions = left_particles.getVariableDataByName<Vecd>("Position");

    BaseParticles &mono_particles = mono_body.body.getBaseParticles();
    const AphiLeftBodyZonePartitionMetrics zone_metrics = hostLeftBodyZonePartitionMetrics(
        left_particles, mono_particles, names.solution, names.solution, names.r_hat, left_positions,
        left_total_real_particles, body_length, body_height, body_width, core_shell, x_interface,
        interface_band_half_width);
    const AphiLeftBodyGaugeErrorMetrics &gauge_metrics = zone_metrics.full_left_half;
    const AphiLeftBodySplitComparisonMetrics &split_metrics = zone_metrics.core_interior_split;

    const AphiDivAReductionMetrics global_div_a = hostTwoBodyContactDivAReductionMetrics(contact_case, names);
    const AphiDivAReductionMetrics left_div_a = hostTwoBodyContactDivAReductionMetrics(
        contact_case, names, [&](const Vecd &position) {
            return isLeftHalfCoreParticle(position, body_length, body_height, body_width, core_shell, x_interface);
        });

    const bool solver_ok = gmresConvergencePassed(contact_metrics.solver_result, tolerance) &&
                           contact_metrics.global_true_rel < 10.0 * tolerance;
    const bool diagnostic_ok = gauge_metrics.particle_count > 0 && zone_metrics.core_interior.particle_count > 0 &&
                               zone_metrics.interface_band.particle_count > 0 &&
                               split_metrics.matched_mono_particles > 0 && split_metrics.missing_mono_particles == 0;
    const bool passed = solver_ok && diagnostic_ok;

    std::cout << "test_3d_aphi_ck_contact_left_field_error_diagnostic"
              << " contact_global_true_rel=" << contact_metrics.global_true_rel
              << " contact_max_bodywise_true_rel=" << contact_metrics.max_bodywise_true_rel
              << " contact_left_true_rel=" << contact_metrics.left_true_rel
              << " contact_right_true_rel=" << contact_metrics.right_true_rel
              << " contact_left_continuous=" << contact_metrics.left_continuous_error
              << " left_block_L2_rel=" << gauge_metrics.block_L2_rel
              << " left_block_Linf_rel=" << gauge_metrics.block_Linf_rel
              << " left_A_real_mean_error=" << gauge_metrics.A_real_mean_error
              << " left_phi_real_mean_error=" << gauge_metrics.phi_real_mean_error
              << " left_block_L2_after_mean_sub=" << gauge_metrics.block_L2_after_mean_subtraction_rel
              << " left_block_Linf_after_mean_sub=" << gauge_metrics.block_Linf_after_mean_subtraction_rel
              << " core_interior_count=" << zone_metrics.core_interior.particle_count
              << " core_interior_L2_rel=" << zone_metrics.core_interior.block_L2_rel
              << " core_interior_Linf_rel=" << zone_metrics.core_interior.block_Linf_rel
              << " core_interior_L2_after_mean_sub=" << zone_metrics.core_interior.block_L2_after_mean_subtraction_rel
              << " core_interior_contact_vs_exact_L2=" << zone_metrics.core_interior_split.contact_vs_exact_L2_rel
              << " core_interior_contact_vs_mono_L2=" << zone_metrics.core_interior_split.contact_vs_mono_L2_rel
              << " interface_band_count=" << zone_metrics.interface_band.particle_count
              << " interface_band_L2_rel=" << zone_metrics.interface_band.block_L2_rel
              << " interface_band_Linf_rel=" << zone_metrics.interface_band.block_Linf_rel
              << " interface_band_L2_after_mean_sub="
              << zone_metrics.interface_band.block_L2_after_mean_subtraction_rel
              << " interface_band_contact_vs_exact_L2=" << zone_metrics.interface_band_split.contact_vs_exact_L2_rel
              << " interface_band_contact_vs_mono_L2=" << zone_metrics.interface_band_split.contact_vs_mono_L2_rel
              << " div_A_relative=" << global_div_a.div_a_relative
              << " div_A_Linf=" << global_div_a.div_a_Linf << " left_div_A_relative=" << left_div_a.div_a_relative
              << " div_A_level=" << divAGaugeDiagnosticLevel(global_div_a.div_a_relative)
              << " contact_vs_exact_L2_rel=" << split_metrics.contact_vs_exact_L2_rel
              << " contact_vs_exact_Linf_rel=" << split_metrics.contact_vs_exact_Linf_rel
              << " contact_vs_mono_L2_rel=" << split_metrics.contact_vs_mono_L2_rel
              << " contact_vs_mono_Linf_rel=" << split_metrics.contact_vs_mono_Linf_rel
              << " mono_vs_exact_L2_rel=" << split_metrics.mono_vs_exact_L2_rel
              << " mono_vs_exact_Linf_rel=" << split_metrics.mono_vs_exact_Linf_rel
              << " matched_mono_particles=" << split_metrics.matched_mono_particles
              << " missing_mono_particles=" << split_metrics.missing_mono_particles
              << " left_particle_count=" << gauge_metrics.particle_count << " solver_ok=" << (solver_ok ? 1 : 0)
              << " diagnostic_ok=" << (diagnostic_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
