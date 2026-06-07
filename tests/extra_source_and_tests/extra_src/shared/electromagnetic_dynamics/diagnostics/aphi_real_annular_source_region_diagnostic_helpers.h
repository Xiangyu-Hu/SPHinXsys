#ifndef APHI_REAL_ANNULAR_SOURCE_REGION_DIAGNOSTIC_HELPERS_H
#define APHI_REAL_ANNULAR_SOURCE_REGION_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/benchmark/aphi_team7_canonical_case_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h"

#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiRealAnnularSourceRegionSpec
{
    Real dp = 0.1;
    Real body_length = benchmark::AphiTeam7PhysicalDimensions::length;
    Real body_height = benchmark::AphiTeam7PhysicalDimensions::height;
    Real body_width = benchmark::AphiTeam7PhysicalDimensions::width;
    benchmark::AphiAnnularSourceGeometry annulus{};
    Real omega = benchmark::AphiTeam7CanonicalCaseSpec::omega;
    Real phi_gauge_penalty = benchmark::AphiTeam7CanonicalCaseSpec::phi_gauge_penalty;
    Real impressed_current_amplitude = benchmark::AphiTeam7CanonicalCaseSpec::impressed_current_amplitude;
    Real tolerance = benchmark::AphiTeam7CanonicalCaseSpec::tolerance;
    UnsignedInt restart_dimension = benchmark::AphiTeam7CanonicalCaseSpec::restart_dimension;
    UnsignedInt max_outer_iterations = benchmark::AphiTeam7CanonicalCaseSpec::max_outer_iterations;
    Real min_conductor_J = 1.0e-12;
    Real min_conductor_joule_integral = 1.0e-14;
    Real source_to_conductor_joule_ratio_cap = 1.0e-6;
    bool write_vtp = false;
};

inline AphiSourceDrivenEmSolveMetrics runRealAnnularSourceRegionDiagnostic(int ac, char *av[],
                                                                           const AphiRealAnnularSourceRegionSpec &spec)
{
    const Real boundary_width = 3.0 * spec.dp;
    const benchmark::AphiTeam7LikeUnitBoxLayout layout =
        benchmark::buildTeam7LayoutForBox(spec.body_length, spec.body_height, spec.body_width);
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    AphiJouleHeatingFieldNames joule_fields;
    AphiSourceDrivenEmSolveFieldNames obs_fields;
    AphiVariableNames names;

    AphiLhsTestBody test_body(spec.dp, spec.body_length, spec.body_height, spec.body_width, boundary_width, ac, av);
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    registerSourceDrivenEmSolveObservableFields(test_body.body, obs_fields);

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, benchmark::AssignTeam7LikeRegionMaterialsCK> assign_material(
        test_body.body, layout, names.material);
    StateDynamics<MainExecutionPolicy, AssignTeam7MaterialRegionIdCK> assign_region_id(
        test_body.body, layout, obs_fields.material_region_id);
    StateDynamics<MainExecutionPolicy, benchmark::AssignZeroSigmaInAnnularRegionCK> assign_annular_sigma_zero(
        test_body.body, spec.annulus, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignImpressedCurrentRhsAnnularCK> assign_annular_rhs(
        test_body.body, names.rhs, spec.annulus, coil_current_real, coil_current_imag, spec.impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = spec.omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = spec.phi_gauge_penalty;
    options.use_a_divergence_penalty = false;

    AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(spec.tolerance, spec.restart_dimension, spec.max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options, solver_options);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, spec.omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);

    (void)register_joule_fields;
    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_region_id.exec();
    assign_annular_sigma_zero.exec();
    assign_annular_rhs.exec();
    test_body.updateRelations();
    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult em_result = solver.solve();

    execBodyCurlBFromADiagnostic(test_body.body, test_body.inner(), names, obs_fields.b_real, obs_fields.b_imag,
                                 AphiBCurlDiagnosticMode::BCorrectedGrad);
    StateDynamics<MainExecutionPolicy, AphiComplexVecdMagnitudeCK> j_magnitude(
        test_body.body, joule_fields.current_density_real, joule_fields.current_density_imag, obs_fields.j_magnitude);
    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    j_magnitude.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    const auto in_annular_source = [&](const Vecd &position) {
        return benchmark::insideAnnularSourceRegion(position, spec.annulus);
    };
    const auto in_physical_conductor = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Conductor);
    };
    const auto in_physical_air = [&](const Vecd &position) {
        return team7ParticleInPhysicalRegion(position, layout, spec.body_length, spec.body_height, spec.body_width,
                                             AphiBenchmarkMaterialRegion::Air);
    };

    AphiSourceDrivenEmSolveMetrics metrics;
    metrics.particles = total_real_particles;
    metrics.num_iterations = em_result.outer_iteration_count;
    metrics.final_residual = em_result.final_true_relative_residual;
    metrics.converged = gmresConvergencePassed(em_result, spec.tolerance);
    metrics.conductor_Joule_integral = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_physical_conductor, joule_fields.joule_heat_source);
    metrics.source_Joule_integral = hostParticleRegionVolWeightedJoulePower(
        particles, positions, total_real_particles, in_annular_source, joule_fields.joule_heat_source);
    metrics.air_Joule_integral = hostParticleRegionVolWeightedJoulePower(particles, positions, total_real_particles,
                                                                         in_physical_air, joule_fields.joule_heat_source);
    metrics.source_rhs_l2_source_region =
        hostParticleRegionBlockNorm(particles, names.rhs, positions, total_real_particles, in_annular_source);
    metrics.source_rhs_l2 = hostBlockNorm(particles, names.rhs, total_real_particles);
    metrics.max_abs_J_conductor = hostParticleRegionScalarMax(particles, obs_fields.j_magnitude, positions,
                                                              total_real_particles, in_physical_conductor);
    metrics.finite_field_check =
        hostSourceDrivenFieldsFinite(particles, names, joule_fields, obs_fields, total_real_particles);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (benchmark::insideAnnularSourceRegion(positions[i], spec.annulus))
        {
            metrics.particle_count_source += 1;
        }
        if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Conductor))
        {
            metrics.particle_count_conductor += 1;
        }
        else if (team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Coil))
        {
            (void)0;
        }
        else
        {
            metrics.particle_count_air += 1;
        }
    }

    if (spec.write_vtp)
    {
        writeSourceDrivenEmSolveVtp(test_body.sph_system, test_body.body, names, joule_fields, obs_fields, names.material);
    }
    return metrics;
}

inline bool realAnnularSourceRegionDiagnosticPassed(const AphiSourceDrivenEmSolveMetrics &metrics,
                                                    const AphiRealAnnularSourceRegionSpec &spec)
{
    const Real max_source_joule =
        std::max(Real(1.0e-12), spec.source_to_conductor_joule_ratio_cap * metrics.conductor_Joule_integral);
    return metrics.converged && metrics.finite_field_check && metrics.particle_count_source > 0 &&
           metrics.particle_count_conductor > 0 && metrics.source_rhs_l2_source_region > 0.0 &&
           metrics.source_Joule_integral <= max_source_joule &&
           metrics.max_abs_J_conductor > spec.min_conductor_J &&
           metrics.conductor_Joule_integral > spec.min_conductor_joule_integral;
}

inline void printRealAnnularSourceRegionDiagnosticMetrics(const char *test_name,
                                                          const AphiSourceDrivenEmSolveMetrics &metrics,
                                                          const AphiRealAnnularSourceRegionSpec &spec, bool passed)
{
    printSourceDrivenEmSolveMetrics(test_name, metrics, passed);
    const Real source_to_conductor_joule =
        metrics.source_Joule_integral / (metrics.conductor_Joule_integral + TinyReal);
    std::cout << test_name << " source_region_model=real_annular annular_geometry_implemented=1"
              << " annulus_center=(" << spec.annulus.center_x << "," << spec.annulus.center_y << ","
              << spec.annulus.center_z << ")"
              << " annulus_inner_radius=" << spec.annulus.inner_radius
              << " annulus_outer_radius=" << spec.annulus.outer_radius
              << " annulus_z_half_height=" << spec.annulus.z_half_height
              << " source_to_conductor_joule=" << source_to_conductor_joule
              << " source_rhs_l2_source_region=" << metrics.source_rhs_l2_source_region << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif
