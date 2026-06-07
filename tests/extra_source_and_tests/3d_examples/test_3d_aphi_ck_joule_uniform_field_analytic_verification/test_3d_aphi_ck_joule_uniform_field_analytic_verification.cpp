/**
 * Stage 10A-1: analytic Joule verification — uniform sigma, A=0, phi=-E0*x.
 * Expected q = 0.5 * sigma * E0^2 (no GMRES; validates grad/E/Joule pipeline only).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"

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
    const Real core_shell = 2.0 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real electric_field_x = 3.0;
    const Real expected_joule = Real(0.5) * sigma * electric_field_x * electric_field_x;
    const Real max_core_mean_relative_error = 0.10;

    const AphiJouleHeatingFieldNames joule_fields;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    RegisterAphiJouleHeatingFieldsCK register_joule_fields(test_body.body, joule_fields);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu,
                                                                                             names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu,
                                                                                 names.material);
    StateDynamics<MainExecutionPolicy, AssignUniformLinearPhiFieldsCK> assign_uniform_phi(
        test_body.body, electric_field_x, names.solution);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, omega, names.solution, joule_fields);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule_source(
        test_body.body, names.material, joule_fields);

    (void)register_joule_fields;
    initialize_aphi_variables.exec();
    set_material.exec();
    assign_uniform_phi.exec();
    test_body.updateRelations();

    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule_source.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Real core_mean_joule =
        hostCoreMeanScalar(particles, joule_fields.joule_heat_source, total_real_particles, body_length, body_height,
                           body_width, core_shell);
    const Real core_mean_relative_error = relativeMetricChange(expected_joule, core_mean_joule);

    syncVariableToHost<Vecd>(particles, joule_fields.electric_field_a_real);
    const Vecd *e_real = particles.getVariableDataByName<Vecd>(joule_fields.electric_field_a_real);
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    Real core_mean_ex = 0.0;
    size_t core_count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        core_mean_ex += e_real[i][0];
        ++core_count;
    }
    core_mean_ex /= static_cast<Real>(core_count) + TinyReal;

    const bool passed = core_mean_relative_error <= max_core_mean_relative_error;

    std::cout << "test_3d_aphi_ck_joule_uniform_field_analytic_verification"
              << " particles=" << total_real_particles << " E0=" << electric_field_x
              << " expected_joule=" << expected_joule << " core_mean_joule=" << core_mean_joule
              << " core_mean_rel_error=" << core_mean_relative_error << " core_mean_Ex=" << core_mean_ex
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
