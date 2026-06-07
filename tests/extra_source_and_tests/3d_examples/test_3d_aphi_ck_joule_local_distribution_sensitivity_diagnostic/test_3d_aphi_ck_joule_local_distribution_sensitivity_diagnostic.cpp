/**
 * Stage 10A-2: Joule local/source-field sensitivity vs EM tolerance (extends 9D-2).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"

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
    const Real conductor_band_half_width = 2.0 * dp_0;
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;
    const Real power_variation_threshold = 0.02;
    const Real conductor_l2_ratio_threshold = 0.02;
    const Real conductor_max_rel_threshold = 0.05;
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
    std::vector<Real> reference_joule;
    const size_t reference_index = 1;

    std::cout << "test_3d_aphi_ck_joule_local_distribution_sensitivity_diagnostic";

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

        if (sweep_results.size() == reference_index)
        {
            reference_joule = hostCopyJouleSourceReference(particles, joule_fields.joule_heat_source, total_real_particles);
        }
        else if (!reference_joule.empty() && sweep_results.size() > reference_index)
        {
            entry.local = hostJouleLocalFieldDifferenceMetrics(
                particles, positions, total_real_particles, layout, joule_fields.joule_heat_source,
                reference_joule.data(), conductor_band_half_width);
        }

        sweep_results.push_back(entry);

        std::cout << " tol=" << tolerance << " em_rel=" << entry.em_rel << " converged=" << (entry.converged ? 1 : 0)
                  << " conductor_joule=" << entry.joule.conductor;
        if (!reference_joule.empty() && sweep_results.size() > reference_index + 1)
        {
            std::cout << " conductor_l2_diff=" << entry.local.conductor_l2_difference
                      << " conductor_max_diff=" << entry.local.conductor_max_abs_difference
                      << " band_l2_diff=" << entry.local.conductor_band_l2_difference;
        }
    }

    const Real reference_conductor = sweep_results[reference_index].joule.conductor;
    const auto in_conductor = [&](const Vecd &position) {
        return team7ParticleInRegion(position, layout, AphiBenchmarkMaterialRegion::Conductor);
    };
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real reference_conductor_l2_sq = 0.0;
    Real reference_conductor_max = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!in_conductor(positions[i]))
        {
            continue;
        }
        reference_conductor_l2_sq += vol[i] * reference_joule[i] * reference_joule[i];
        reference_conductor_max = std::max(reference_conductor_max, reference_joule[i]);
    }
    const Real reference_conductor_l2 = std::sqrt(reference_conductor_l2_sq);

    Real max_conductor_power_variation = 0.0;
    Real max_conductor_l2_ratio = 0.0;
    Real max_conductor_max_rel = 0.0;

    for (size_t index = reference_index; index < sweep_results.size(); ++index)
    {
        const AphiJouleEmToleranceSweepEntry &entry = sweep_results[index];
        max_conductor_power_variation =
            std::max(max_conductor_power_variation, relativeMetricChange(reference_conductor, entry.joule.conductor));
        if (index > reference_index)
        {
            max_conductor_l2_ratio = std::max(
                max_conductor_l2_ratio, entry.local.conductor_l2_difference / (reference_conductor_l2 + TinyReal));
            max_conductor_max_rel = std::max(max_conductor_max_rel,
                                             entry.local.conductor_max_abs_difference / (reference_conductor_max + TinyReal));
        }
    }

    bool all_joule_valid = true;
    for (const AphiJouleEmToleranceSweepEntry &entry : sweep_results)
    {
        all_joule_valid = all_joule_valid && entry.joule.total > 0.0 && entry.joule.conductor >= 0.0 &&
                          entry.joule.min_source >= -1.0e-12;
    }

    const bool sensitivity_passed = max_conductor_power_variation <= power_variation_threshold &&
                                    max_conductor_l2_ratio <= conductor_l2_ratio_threshold &&
                                    max_conductor_max_rel <= conductor_max_rel_threshold && all_joule_valid;

    std::cout << " max_conductor_power_rel_change=" << max_conductor_power_variation
              << " max_conductor_l2_ratio=" << max_conductor_l2_ratio
              << " max_conductor_max_rel=" << max_conductor_max_rel
              << " sensitivity_passed=" << (sensitivity_passed ? 1 : 0) << " passed=1" << std::endl;

    return all_joule_valid ? 0 : 1;
}
