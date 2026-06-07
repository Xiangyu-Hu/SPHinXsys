#ifndef APHI_CONTACT_SOURCE_DRIVEN_HEATING_BASELINE_HELPERS_H
#define APHI_CONTACT_SOURCE_DRIVEN_HEATING_BASELINE_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_contact_a_divergence_penalty_two_body_mms_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_physical_region_audit_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.hpp"
#include "electromagnetic_dynamics/aphi_multibody_contact_gmres_ck.h"

#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactSourceDrivenHeatingBaselineMetrics
{
    bool converged = false;
    UnsignedInt num_iterations = 0;
    Real final_residual = 0.0;
    Real plate_joule_integral = 0.0;
    Real air_joule_integral = 0.0;
    Real max_abs_E = 0.0;
    Real max_abs_J = 0.0;
    bool finite_fields = true;
    AphiConductorInterfaceSpikeMetrics interface_spike{};
    AphiThreeBodyContactAuditSummary three_body_audit{};
    AphiBodyRegionAuditMetrics left_audit{};
    AphiBodyRegionAuditMetrics right_audit{};
    Real total_source_rhs_l2 = 0.0;
    Real total_conductor_joule = 0.0;
    Real total_air_joule = 0.0;
    bool is_three_body = false;
};

inline void initializeTwoBodySourceDrivenContactFields(AphiTwoBodyInterfaceCase &case_setup,
                                                       const benchmark::AphiTeam7LikeUnitBoxLayout &layout,
                                                       const AphiVariableNames &names)
{
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(
        case_setup.left_body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(
        case_setup.right_body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, benchmark::AssignTeam7LikeRegionMaterialsCK> assign_left_materials(
        case_setup.left_body, layout, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignTeam7LikeRegionMaterialsCK> assign_right_materials(
        case_setup.right_body, layout, names.material);

    initialize_left.exec();
    initialize_right.exec();
    assign_left_materials.exec();
    assign_right_materials.exec();
}

/** P4: air (left) + TEAM7 conductor/coil/air (right), coupled Contact GMRES, lambda_A off. */
inline AphiContactSourceDrivenHeatingBaselineMetrics runTwoBodyContactSourceDrivenHeatingBaseline(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = benchmark::AphiTeam7PhysicalDimensions::length;
    const Real body_height = benchmark::AphiTeam7PhysicalDimensions::height;
    const Real body_width = benchmark::AphiTeam7PhysicalDimensions::width;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real omega = benchmark::AphiTeam7CanonicalCaseSpec::omega;
    const Real phi_gauge_penalty = benchmark::AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    const Real tolerance = benchmark::AphiTeam7CanonicalCaseSpec::tolerance;
    const Real impressed_current_amplitude = benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    const benchmark::AphiTeam7LikeUnitBoxLayout layout =
        benchmark::buildTeam7LayoutForBox(body_length, body_height, body_width);

    AphiTwoBodyInterfaceCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = false;
    options.a_divergence_penalty = 0.0;

    initializeTwoBodySourceDrivenContactFields(case_setup, layout, names);

    const AphiMatrixFreeSolverOptions solver_options =
        defaultCoupledContactGMRESOptions(tolerance, benchmark::AphiTeam7CanonicalCaseSpec::restart_dimension,
                                        benchmark::AphiTeam7CanonicalCaseSpec::max_outer_iterations);
    const AphiGMRESWorkspaceNames gmres_workspace =
        buildAphiGMRESWorkspaceNames(solver_options.gmres.restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_left_workspace(case_setup.left_body, solver_options.gmres.restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_right_workspace(case_setup.right_body, solver_options.gmres.restart_dimension);
    RegisterAphiJouleHeatingFieldsCK register_left_joule(case_setup.left_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_right_joule(case_setup.right_body, joule_names);
    (void)register_left_workspace;
    (void)register_right_workspace;
    (void)register_left_joule;
    (void)register_right_joule;

    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left_rhs(case_setup.left_body, names.rhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right_rhs(case_setup.right_body, names.rhs);
    StateDynamics<MainExecutionPolicy, benchmark::AssignImpressedCurrentRhsCK> assign_coil_source(
        case_setup.left_body, names.rhs, layout.coil, coil_current_real, coil_current_imag, impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left_solution(case_setup.left_body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right_solution(case_setup.right_body, names.solution);

    zero_left_rhs.exec();
    zero_right_rhs.exec();
    assign_coil_source.exec();
    zero_left_solution.exec();
    zero_right_solution.exec();
    case_setup.updateRelations();

    const StdVec<AphiMultiBodyContactEntry> entries = buildTwoBodyMultiBodyEntries(case_setup);
    AphiMultiBodyContactGMRESSolverCK<MainExecutionPolicy> solver(entries, names, gmres_workspace, options,
                                                                  solver_options.gmres);
    const AphiMatrixFreeSolverResult solver_result = toMatrixFreeSolverResult(solver.solve());
    case_setup.updateRelations();

    execTwoBodyContactJoulePostProcess(case_setup, names, omega, joule_names);

    const AphiTeam7PlateObservables plate_observables = hostTeam7PlateObservablesFromBody(
        case_setup.right_body, names, joule_names, layout, body_length, body_height, body_width, core_shell);

    BaseParticles &left_particles = case_setup.left_body.getBaseParticles();
    const size_t left_particles_count = left_particles.TotalRealParticles();
    syncVariableToHost<Vecd>(left_particles, "Position");
    const Vecd *left_positions = left_particles.getVariableDataByName<Vecd>("Position");
    const auto all_left_particles = [](const Vecd &) { return true; };
    const Real left_joule_integral = hostParticleRegionVolWeightedJoulePower(
        left_particles, left_positions, left_particles_count, all_left_particles, joule_names.joule_heat_source);

    AphiContactSourceDrivenHeatingBaselineMetrics metrics;
    metrics.converged = gmresConvergencePassed(solver_result, tolerance);
    metrics.num_iterations = solver_result.outer_iteration_count;
    metrics.final_residual = solver_result.final_true_relative_residual;
    metrics.plate_joule_integral = plate_observables.plate_joule_power;
    metrics.air_joule_integral = left_joule_integral;
    metrics.max_abs_E = plate_observables.plate_E_L2;
    metrics.max_abs_J = plate_observables.plate_j_L2;
    metrics.finite_fields = std::isfinite(metrics.plate_joule_integral) && metrics.plate_joule_integral > 0.0 &&
                          std::isfinite(metrics.air_joule_integral);

    syncVariableToHost<Vecd>(case_setup.right_body.getBaseParticles(), "Position");
    const Vecd *right_positions = case_setup.right_body.getBaseParticles().getVariableDataByName<Vecd>("Position");
    const size_t right_count = case_setup.right_body.getBaseParticles().TotalRealParticles();
    const Real x_interface = 0.5 * body_length;
    const Real interface_band_half_width = dp_0;
    metrics.interface_spike = hostConductorInterfaceSpikeMetrics(
        case_setup.right_body.getBaseParticles(), joule_names, right_positions, right_count, layout, body_length,
        body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0);
    const Real conductor_sigma = layout.conductor_material.sigma;
    metrics.left_audit = hostAuditBodyRegions(case_setup.left_body, "left", names, joule_names, layout, body_length,
                                              body_height, body_width, conductor_sigma);
    metrics.right_audit = hostAuditBodyRegions(case_setup.right_body, "right", names, joule_names, layout, body_length,
                                               body_height, body_width, conductor_sigma);
    metrics.total_source_rhs_l2 = metrics.left_audit.rhs_l2;
    metrics.total_conductor_joule = metrics.right_audit.joule_integral_conductor;
    metrics.total_air_joule = metrics.left_audit.joule_integral_air + metrics.right_audit.joule_integral_air;
    metrics.plate_joule_integral = metrics.total_conductor_joule;
    metrics.air_joule_integral = metrics.total_air_joule;
    metrics.is_three_body = false;
    return metrics;
}

/** P4b: three-body TEAM7-like Contact source-driven baseline (lambda_A off). */
inline AphiContactSourceDrivenHeatingBaselineMetrics runThreeBodyContactSourceDrivenHeatingBaseline(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = benchmark::AphiTeam7PhysicalDimensions::length;
    const Real body_height = benchmark::AphiTeam7PhysicalDimensions::height;
    const Real body_width = benchmark::AphiTeam7PhysicalDimensions::width;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real omega = benchmark::AphiTeam7CanonicalCaseSpec::omega;
    const Real phi_gauge_penalty = benchmark::AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    const Real tolerance = benchmark::AphiTeam7CanonicalCaseSpec::tolerance;
    const Real impressed_current_amplitude = benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    AphiTeam7ThreeBodyContactCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    options.use_a_divergence_penalty = false;
    options.a_divergence_penalty = 0.0;

    const AphiMatrixFreeSolverOptions solver_options =
        defaultCoupledContactGMRESOptions(tolerance, benchmark::AphiTeam7CanonicalCaseSpec::restart_dimension,
                                        benchmark::AphiTeam7CanonicalCaseSpec::max_outer_iterations);
    const AphiTeam7ContactCoupledGmresResult result = runTeam7ThreeBodyCoupledContactGmres(
        case_setup, names, joule_names, options, solver_options, body_length, body_height, body_width, core_shell,
        coil_current_real, coil_current_imag, impressed_current_amplitude, true);

    AphiContactSourceDrivenHeatingBaselineMetrics metrics;
    metrics.converged = gmresConvergencePassed(result.solver_result, tolerance);
    metrics.num_iterations = result.solver_result.outer_iteration_count;
    metrics.final_residual = result.solver_result.final_true_relative_residual;
    metrics.plate_joule_integral = result.plate.plate_joule_power;
    metrics.max_abs_E = result.plate.plate_E_L2;
    metrics.max_abs_J = result.plate.plate_j_L2;
    metrics.finite_fields = std::isfinite(metrics.plate_joule_integral) && metrics.plate_joule_integral > 0.0;

    syncVariableToHost<Vecd>(case_setup.plate_body.getBaseParticles(), "Position");
    const Vecd *plate_positions = case_setup.plate_body.getBaseParticles().getVariableDataByName<Vecd>("Position");
    const size_t plate_count = case_setup.plate_body.getBaseParticles().TotalRealParticles();
    const Real x_interface = case_setup.layout.conductor.xmin;
    const Real interface_band_half_width = dp_0;
    metrics.interface_spike = hostConductorInterfaceSpikeMetrics(
        case_setup.plate_body.getBaseParticles(), joule_names, plate_positions, plate_count, case_setup.layout,
        body_length, body_height, body_width, core_shell, x_interface, interface_band_half_width, dp_0);
    metrics.is_three_body = true;
    metrics.three_body_audit =
        hostThreeBodyContactAuditSummary(case_setup, names, joule_names, body_length, body_height, body_width);
    metrics.air_joule_integral = metrics.three_body_audit.global_air_joule;
    return metrics;
}

inline bool contactSourceDrivenHeatingBaselinePassed(const AphiContactSourceDrivenHeatingBaselineMetrics &metrics)
{
    if (!metrics.converged || !metrics.finite_fields || metrics.plate_joule_integral <= 1.0e-6 ||
        metrics.max_abs_J <= 1.0e-8)
    {
        return false;
    }
    if (!metrics.interface_spike.has_bulk_reference || !metrics.interface_spike.spike_hard_pass)
    {
        return false;
    }
    if (!metrics.is_three_body)
    {
        return metrics.total_source_rhs_l2 > 0.0 && metrics.total_conductor_joule > 1.0e-6 &&
               twoBodyContactAuditPassed(metrics.left_audit, metrics.right_audit, 0.05);
    }
    return threeBodyContactAuditPassed(metrics.three_body_audit, 0.05);
}

inline bool contactThreeBodyCompositionPassed(const AphiThreeBodyContactAuditSummary &audit)
{
    return audit.total_source_rhs_l2 > 0.0 && audit.plate.particle_count_conductor > 0 &&
           audit.coil.particle_count_source > 0 && audit.air.particle_count_air > 0 &&
           audit.global_conductor_joule > 1.0e-6 &&
           audit.global_air_joule / (audit.global_conductor_joule + TinyReal) < 0.05 &&
           audit.coil.sigma_max <= 1.0e-14;
}

inline void printContactSourceDrivenHeatingBaselineMetrics(const char *test_name,
                                                           const AphiContactSourceDrivenHeatingBaselineMetrics &metrics,
                                                           bool passed)
{
    std::cout << test_name << " passed=" << (passed ? 1 : 0) << " converged=" << (metrics.converged ? 1 : 0)
              << " num_iterations=" << metrics.num_iterations << " final_residual=" << metrics.final_residual
              << " plate_Joule_integral=" << metrics.plate_joule_integral
              << " air_Joule_integral=" << metrics.air_joule_integral << " max_abs_E=" << metrics.max_abs_E
              << " max_abs_J=" << metrics.max_abs_J << " lambda_A=off" << std::endl;
    printConductorInterfaceSpikeMetrics(test_name, metrics.interface_spike, 10.0, 50.0);
    if (metrics.is_three_body)
    {
        printThreeBodyContactAuditSummary(test_name, metrics.three_body_audit);
        std::cout << test_name << " three_body_audit_passed="
                  << (threeBodyContactAuditPassed(metrics.three_body_audit) ? 1 : 0) << std::endl;
    }
    else
    {
        printBodyRegionAuditMetrics(test_name, metrics.left_audit);
        printBodyRegionAuditMetrics(test_name, metrics.right_audit);
        std::cout << test_name << " total_source_rhs_l2=" << metrics.total_source_rhs_l2
                  << " total_conductor_joule=" << metrics.total_conductor_joule
                  << " total_air_joule=" << metrics.total_air_joule << " two_body_audit_passed="
                  << (twoBodyContactAuditPassed(metrics.left_audit, metrics.right_audit) ? 1 : 0) << std::endl;
    }
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
