#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics;

namespace
{

using MainExecutionPolicy = execution::MainExecutionPolicy;

class VariableRegistrationBoxShape : public ComplexShape
{
  public:
    VariableRegistrationBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

bool vectorFieldIsZero(const Vecd *field, size_t total_real_particles)
{
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (field[i].norm() > 1.0e-12)
        {
            return false;
        }
    }
    return true;
}

bool scalarFieldIsZero(const Real *field, size_t total_real_particles)
{
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (std::abs(field[i]) > 1.0e-12)
        {
            return false;
        }
    }
    return true;
}

bool scalarFieldEquals(const Real *field, size_t total_real_particles, Real reference)
{
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (std::abs(field[i] - reference) > 1.0e-12)
        {
            return false;
        }
    }
    return true;
}

bool blockFieldIsZero(BaseParticles &particles, const AphiBlockNames &block_names, size_t total_real_particles)
{
    return vectorFieldIsZero(particles.getVariableDataByName<Vecd>(block_names.a_real), total_real_particles) &&
           vectorFieldIsZero(particles.getVariableDataByName<Vecd>(block_names.a_imag), total_real_particles) &&
           scalarFieldIsZero(particles.getVariableDataByName<Real>(block_names.phi_real), total_real_particles) &&
           scalarFieldIsZero(particles.getVariableDataByName<Real>(block_names.phi_imag), total_real_particles);
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.2;
    const Vecd halfsize(0.5, 0.5, 0.5);
    const Vecd center(0.5, 0.5, 0.5);
    const Real boundary_width = 2.0 * dp_0;
    const Real default_sigma = 3.5;
    const Real default_nu = 7.25;
    const Real reset_sigma = 2.0;
    const Real reset_nu = 4.0;

    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(1.0 + boundary_width, 1.0 + boundary_width, 1.0 + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    SolidBody body(sph_system, makeShared<VariableRegistrationBoxShape>("AphiVariableBody", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.generateParticles<BaseParticles, Lattice>();

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(body, default_sigma, default_nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> reset_material_properties(body, reset_sigma, reset_nu, names.material);
    initialize_aphi_variables.exec();

    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();

    const Real *sigma = particles.getVariableDataByName<Real>(names.material.sigma);
    const Real *nu = particles.getVariableDataByName<Real>(names.material.nu);

    const bool initialization_passed = blockFieldIsZero(particles, names.solution, total_real_particles) &&
                                       blockFieldIsZero(particles, names.rhs, total_real_particles) &&
                                       blockFieldIsZero(particles, names.lhs, total_real_particles) &&
                                       blockFieldIsZero(particles, names.residual, total_real_particles) &&
                                       blockFieldIsZero(particles, names.r_hat, total_real_particles) &&
                                       blockFieldIsZero(particles, names.search, total_real_particles) &&
                                       blockFieldIsZero(particles, names.v, total_real_particles) &&
                                       blockFieldIsZero(particles, names.s, total_real_particles) &&
                                       blockFieldIsZero(particles, names.t, total_real_particles) &&
                                       scalarFieldEquals(sigma, total_real_particles, default_sigma) &&
                                       scalarFieldEquals(nu, total_real_particles, default_nu);

    reset_material_properties.exec();

    const bool material_reset_passed = scalarFieldEquals(sigma, total_real_particles, reset_sigma) &&
                                       scalarFieldEquals(nu, total_real_particles, reset_nu);
    const bool passed = initialization_passed && material_reset_passed;

    std::cout << "test_3d_aphi_ck_variable_registration"
              << " total_real_particles=" << total_real_particles
              << " sigma=" << sigma[0]
              << " nu=" << nu[0]
              << " material_reset_passed=" << (material_reset_passed ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
