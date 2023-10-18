/**
 * @file 	channel_flow_shell.cpp
 * @brief 	This is a test of fluid-shell interaction.
 * @author 	Weiyi Kong
 */
#include "sphinxsys.h"

#include "case.h" //	case file to setup the test case
using namespace SPH;

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
    //----------------------------------------------------------------------
    // Contact
    //----------------------------------------------------------------------
    InnerRelation shell_inner_contact(shell);
    SimpleDynamics<ShellInitialCurvature> shell_initial_curvature(shell_inner_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    // calculate initial curvature
    shell_initial_curvature.exec();
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    shell.addBodyStateForRecording<Real>("MeanCurvature");
    shell.addBodyStateForRecording<Real>("GaussianCurvature");
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    write_real_body_states.writeToFile();

    return 0;
}
