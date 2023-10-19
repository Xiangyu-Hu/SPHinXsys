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
    const Real end_time = 10.0;
    const Real target_radius_increase = 0.5 * radius;
    const Real moving_v = target_radius_increase / end_time;
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
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson); // dummy material parameters
    shell.generateParticles<ShellParticleGenerator>();
    //----------------------------------------------------------------------
    // Contact
    //----------------------------------------------------------------------
    InnerRelation shell_inner_contact(shell);
    //----------------------------------------------------------------------
    // Solid Algorithms
    //----------------------------------------------------------------------
    /** Corrected configuration. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration>
        corrected_configuration(shell_inner_contact);
    /** Time step size calculation. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> computing_time_step_size(shell);
    /** stress relaxation. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> stress_relaxation_first_half(shell_inner_contact);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> stress_relaxation_second_half(shell_inner_contact);
    /** Control the displacement. */
    SimpleDynamics<ControlDisplacement> dis_control(shell, moving_v);
    /** Constrain the Boundary. */
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vecd>>>
        cylinder_position_damping(0.2, shell_inner_contact, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vecd>>>
        cylinder_rotation_damping(0.2, shell_inner_contact, "AngularVelocity", physical_viscosity);
    /** Compute curvature. */
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> shell_initial_curvature(shell_inner_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();
    // compute initial curvature after B_ is computed
    shell_initial_curvature.compute_initial_curvature();
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    shell.addBodyStateForRecording<Real>("MeanCurvature");
    shell.addBodyStateForRecording<Real>("GaussianCurvature");
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    write_real_body_states.writeToFile();
    /**
     * From here the time stepping begins.
     * Set the starting time.
     */
    GlobalStaticVariables::physical_time_ = 0.0;

    /** Setup time stepping control parameters. */
    int ite = 0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    /** Statistics for computing time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /**
     * Main loop
     */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integral_time = 0.0;
        while (integral_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            stress_relaxation_first_half.exec(dt);
            dis_control.exec(dt);
            cylinder_position_damping.exec(dt);
            cylinder_rotation_damping.exec(dt);
            dis_control.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integral_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }
        TickCount t2 = TickCount::now();
        shell_initial_curvature.exec();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
