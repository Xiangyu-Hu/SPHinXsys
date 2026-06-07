/**
 * Stage 9D-lite: Joule heating data plumbing after TEAM7-like EM solve (informational).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const Real tolerance = 5.0e-4;
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;

    const AphiTeam7LikeUnitBoxLayout layout;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);
    const AphiJouleHeatingFieldNames joule_fields;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, AssignTeam7LikeRegionMaterialsCK> assign_material(test_body.body, layout,
                                                                                         names.material);
    StateDynamics<MainExecutionPolicy, AssignImpressedCurrentRhsCK> assign_coil_source(
        test_body.body, names.rhs, layout.coil, coil_current_real, coil_current_imag, impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                      solver_options);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);

    (void)register_joule_fields;
    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();

    zero_solution.exec();
    test_body.updateRelations();

    const AphiMatrixFreeSolverResult em_result = solver.solve();

    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    const Real total_joule_power =
        hostVolWeightedSum(particles, joule_fields.joule_heat_source, total_real_particles, [](const Vecd &) {
            return true;
        });
    const Real conductor_joule_power =
        hostTeam7RegionJoulePower(particles, positions, total_real_particles, layout,
                                  AphiBenchmarkMaterialRegion::Conductor, joule_fields.joule_heat_source);
    const Real coil_joule_power = hostTeam7RegionJoulePower(particles, positions, total_real_particles, layout,
                                                            AphiBenchmarkMaterialRegion::Coil,
                                                            joule_fields.joule_heat_source);

    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    const Real *joule_source = particles.getVariableDataByName<Real>(joule_fields.joule_heat_source);
    Real min_joule = joule_source[0];
    for (size_t i = 1; i != total_real_particles; ++i)
    {
        min_joule = std::min(min_joule, joule_source[i]);
    }

    const bool passed = em_result.final_true_relative_residual <= tolerance * 10.0 && total_joule_power > 0.0 &&
                        conductor_joule_power >= 0.0 && min_joule >= -1.0e-12;

    std::cout << "test_3d_aphi_ck_joule_heating_plumbing_diagnostic"
              << " em_rel=" << em_result.final_true_relative_residual << " total_joule_power=" << total_joule_power
              << " conductor_joule_power=" << conductor_joule_power << " coil_joule_power=" << coil_joule_power
              << " min_joule_source=" << min_joule << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
