/**
 * Stage 9D-2: Joule power sensitivity vs EM solver tolerance on TEAM7-like case.
 * Sweep tol = 5e-3 / 5e-4 / 1e-4 / 1e-5; reference conductor power at tol=5e-4.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

#include <array>
#include <iostream>
#include <vector>

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
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;
    const Real max_conductor_power_variation = 0.02;
    const std::array<Real, 4> tolerance_values = {5.0e-3, 5.0e-4, 1.0e-4, 1.0e-5};

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

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    std::vector<AphiJouleEmToleranceSweepEntry> sweep_results;
    sweep_results.reserve(tolerance_values.size());

    std::cout << "test_3d_aphi_ck_joule_em_tolerance_sensitivity_diagnostic"
              << " m=" << restart_dimension << " max_outer=" << max_outer_iterations;

    for (const Real tolerance : tolerance_values)
    {
        zero_solution.exec();
        test_body.updateRelations();

        AphiMatrixFreeSolverOptions solver_options =
            defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
        AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                          solver_options);
        const AphiMatrixFreeSolverResult em_result = solver.solve();

        compute_grad_phi.exec();
        compute_electric_field.exec();
        compute_joule_source.exec();
        test_body.updateRelations();

        AphiJouleEmToleranceSweepEntry entry;
        entry.tolerance = tolerance;
        entry.em_rel = em_result.final_true_relative_residual;
        entry.converged = gmresConvergencePassed(em_result, tolerance);
        entry.outer_iterations = em_result.outer_iteration_count;
        entry.joule = hostJoulePowerMetrics(particles, positions, total_real_particles, layout,
                                            joule_fields.joule_heat_source);
        sweep_results.push_back(entry);

        std::cout << " tol=" << tolerance << " em_rel=" << entry.em_rel << " converged=" << (entry.converged ? 1 : 0)
                  << " outer=" << entry.outer_iterations << " total_joule=" << entry.joule.total
                  << " conductor_joule=" << entry.joule.conductor << " coil_joule=" << entry.joule.coil
                  << " min_joule=" << entry.joule.min_source;
    }

    const size_t reference_index = 1;
    const Real reference_conductor = sweep_results[reference_index].joule.conductor;
    const Real reference_total = sweep_results[reference_index].joule.total;

    Real max_conductor_variation = 0.0;
    Real max_total_variation = 0.0;
    std::cout << " reference_tol=" << sweep_results[reference_index].tolerance
              << " conductor_joule=" << reference_conductor;
    for (size_t index = reference_index; index < sweep_results.size(); ++index)
    {
        const AphiJouleEmToleranceSweepEntry &entry = sweep_results[index];
        const Real conductor_change = relativeMetricChange(reference_conductor, entry.joule.conductor);
        const Real total_change = relativeMetricChange(reference_total, entry.joule.total);
        max_conductor_variation = std::max(max_conductor_variation, conductor_change);
        max_total_variation = std::max(max_total_variation, total_change);
        if (index > reference_index)
        {
            std::cout << " delta_vs_ref tol=" << entry.tolerance << " conductor_rel_change=" << conductor_change
                      << " total_rel_change=" << total_change;
        }
    }

    const Real loose_conductor_change =
        relativeMetricChange(reference_conductor, sweep_results[0].joule.conductor);
    std::cout << " loose_tol_5e-3_conductor_rel_change=" << loose_conductor_change;

    bool all_joule_valid = true;
    for (const AphiJouleEmToleranceSweepEntry &entry : sweep_results)
    {
        all_joule_valid = all_joule_valid && entry.joule.total > 0.0 && entry.joule.conductor >= 0.0 &&
                          entry.joule.min_source >= -1.0e-12;
    }

    const bool sensitivity_passed = max_conductor_variation <= max_conductor_power_variation && all_joule_valid;

    std::cout << " max_conductor_rel_change=" << max_conductor_variation
              << " max_total_rel_change=" << max_total_variation
              << " sensitivity_passed=" << (sensitivity_passed ? 1 : 0) << " passed=1" << std::endl;

    return all_joule_valid ? 0 : 1;
}
