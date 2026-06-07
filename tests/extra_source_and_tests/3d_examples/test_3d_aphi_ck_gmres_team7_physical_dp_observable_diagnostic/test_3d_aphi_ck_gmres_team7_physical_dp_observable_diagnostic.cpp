/**
 * Stage 10A-2: physical TEAM7-dimension dp sweep with observable metrics (informational).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"

#include <algorithm>
#include <array>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

struct PhysicalDpObservableResult
{
    Real dp = 0.0;
    size_t particles = 0;
    Real em_rel = 0.0;
    UnsignedInt outer = 0;
    bool converged = false;
    Real conductor_solution_block_max = 0.0;
    Real conductor_joule_power = 0.0;
    Real conductor_joule_max = 0.0;
};

PhysicalDpObservableResult runPhysicalDpCase(int ac, char *av[], Real dp_0)
{
    const Real body_length = AphiTeam7PhysicalDimensions::length;
    const Real body_height = AphiTeam7PhysicalDimensions::height;
    const Real body_width = AphiTeam7PhysicalDimensions::width;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell =
        std::min(Real(2.0) * dp_0, Real(0.125) * std::min({body_length, body_height, body_width}));
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const Real tolerance = 5.0e-4;
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;

    const AphiTeam7LikeUnitBoxLayout layout = buildTeam7LayoutForBox(body_length, body_height, body_width);
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);
    const AphiJouleHeatingFieldNames joule_fields;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);

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

    PhysicalDpObservableResult result;
    result.dp = dp_0;
    result.particles = total_real_particles;
    result.em_rel = em_result.final_true_relative_residual;
    result.outer = em_result.outer_iteration_count;
    result.converged = gmresConvergencePassed(em_result, tolerance);
    result.conductor_solution_block_max =
        coreBlockMaxAbsNormInRegion(particles, names.solution, positions, total_real_particles, body_length,
                                    body_height, body_width, core_shell, layout.conductor);
    result.conductor_joule_power =
        hostTeam7RegionJoulePower(particles, positions, total_real_particles, layout,
                                  AphiBenchmarkMaterialRegion::Conductor, joule_fields.joule_heat_source);

    syncVariableToHost<Real>(particles, joule_fields.joule_heat_source);
    const Real *joule_source = particles.getVariableDataByName<Real>(joule_fields.joule_heat_source);
    result.conductor_joule_max = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!team7ParticleInRegion(positions[i], layout, AphiBenchmarkMaterialRegion::Conductor))
        {
            continue;
        }
        result.conductor_joule_max = std::max(result.conductor_joule_max, joule_source[i]);
    }
    return result;
}

void printResult(const PhysicalDpObservableResult &result)
{
    std::cout << " dp=" << result.dp << " particles=" << result.particles << " em_rel=" << result.em_rel
              << " outer=" << result.outer << " converged=" << (result.converged ? 1 : 0)
              << " conductor_solution_block_max=" << result.conductor_solution_block_max
              << " conductor_joule_power=" << result.conductor_joule_power
              << " conductor_joule_max=" << result.conductor_joule_max;
}

} // namespace

int main(int ac, char *av[])
{
    const std::array<Real, 3> dp_values = {0.15, 0.1, 0.075};

    std::cout << "test_3d_aphi_ck_gmres_team7_physical_dp_observable_diagnostic";
    for (const Real dp : dp_values)
    {
        printResult(runPhysicalDpCase(ac, av, dp));
    }
    std::cout << " passed=1" << std::endl;
    return 0;
}
