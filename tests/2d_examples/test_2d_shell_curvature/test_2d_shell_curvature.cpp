/**
 * @file 	test_2d_shell_curvature.cpp
 * @brief 	This is the case study the curvature of shell.
 * @author  Weiyi Kong
 */

#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real radius_x = 2;
const Real radius_y = 1;
const int particle_number_mid_surface = 100;
const Real resolution_ref = M_PI * (radius_x + radius_y) / Real(particle_number_mid_surface);
const Real shell_thickness = resolution_ref;
std::vector<Vec2d> observer_position = {Vec2d(-radius_x, 0.0), Vec2d(0.0, radius_y)};
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-radius_x, -radius_y), Vec2d(radius_x, radius_y));
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** Particle generator and constraint boundary for shell baffle. */
// x=R*cos(theta), y=R*sin(theta)
class ShellParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Real theta = 0;
        while (theta < 2 * Pi)
        {
            Real x = radius_x * cos(theta);
            Real y = radius_y * sin(theta);
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_ref);
            Vec2d normal_direction = Vec2d(x / radius_x / radius_x, y / radius_y / radius_y);
            normal_direction /= normal_direction.norm();
            initializeSurfaceProperties(normal_direction, shell_thickness);
            Real ds = radius_x * radius_y * sqrt(pow(x / radius_x / radius_x, 2) + pow(y / radius_y / radius_y, 2));
            theta += resolution_ref / ds;
        }
    }
};

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //---------------------------------------------------------------------
    SolidBody shell(sph_system, makeShared<DefaultShape>("Shell"));
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0); // dummy material parameters
    shell.generateParticles<ShellParticleGenerator>();

    ObserverBody observer_curvature(sph_system, "observer_curvature");
    observer_curvature.generateParticles<ParticleGeneratorObserver>(observer_position);
    //----------------------------------------------------------------------
    // Contact
    //----------------------------------------------------------------------
    InnerRelation shell_inner_contact(shell);
    ContactRelation observer_contact(observer_curvature, {&shell});
    //----------------------------------------------------------------------
    // Algorithms
    //----------------------------------------------------------------------
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(shell_inner_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
    ObservingAQuantity<Real> observe_curvature(observer_contact, "Average1stPrincipleCurvature");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    // calculate initial curvature
    shell_curvature.exec();
    //----------------------------------------------------------------------
    //	Record
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    observe_curvature.exec();
    //----------------------------------------------------------------------
    //	gtest
    //----------------------------------------------------------------------
    const Real x_curvature_analytical = radius_x / radius_y / radius_y;
    const Real y_curvature_analytical = radius_y / radius_x / radius_x;
    const Real x_curvature = (*observe_curvature.getParticles()->getVariableByName<Real>("Average1stPrincipleCurvature"))[0];
    const Real y_curvature = (*observe_curvature.getParticles()->getVariableByName<Real>("Average1stPrincipleCurvature"))[1];

    std::cout << "Analytical curvature x: " << x_curvature_analytical
              << "\t Curvature x: " << x_curvature
              << "\t Error: " << std::abs((x_curvature_analytical - x_curvature) / x_curvature_analytical) * 100 << "%"
              << std::endl;

    std::cout << "Analytical curvature y: " << y_curvature_analytical
              << "\t Curvature y: " << y_curvature
              << "\t Error: " << std::abs((y_curvature_analytical - y_curvature) / y_curvature_analytical) * 100 << "%"
              << std::endl;

    // gtest
    EXPECT_NEAR(x_curvature_analytical, x_curvature, x_curvature_analytical * 10e-2); // below 10%
    EXPECT_NEAR(y_curvature_analytical, y_curvature, y_curvature_analytical * 10e-2); // below 10%

    return 0;
}