#include "electromagnetic_dynamics/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real validation_threshold = 1.0e-10;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width);
    test_body.sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, 1.0, 1.0, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, 1.0, 1.0, names.material);

    AphiLhsAssemblyOptions options;
    options.omega = 1.0;

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble(test_body.body, test_body.inner_ck, names, options);

    initialize_aphi_variables.exec();
    set_material.exec();
    test_body.updateRelations();
    assemble.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *lhs_a_real = particles.getVariableDataByName<Vecd>(names.lhs.a_real);
    const Vecd *lhs_a_imag = particles.getVariableDataByName<Vecd>(names.lhs.a_imag);
    const Real *lhs_phi_real = particles.getVariableDataByName<Real>(names.lhs.phi_real);
    const Real *lhs_phi_imag = particles.getVariableDataByName<Real>(names.lhs.phi_imag);

    Real core_max_error = 0.0;
    size_t core_particles = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        ++core_particles;
        core_max_error = std::max(core_max_error, lhs_a_real[i].norm());
        core_max_error = std::max(core_max_error, lhs_a_imag[i].norm());
        core_max_error = std::max(core_max_error, std::abs(lhs_phi_real[i]));
        core_max_error = std::max(core_max_error, std::abs(lhs_phi_imag[i]));
    }

    const bool passed = core_particles > 0 && core_max_error < validation_threshold;

    std::cout << "test_3d_aphi_ck_lhs_full_zero"
              << " core_particles=" << core_particles
              << " core_max_error=" << core_max_error
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
