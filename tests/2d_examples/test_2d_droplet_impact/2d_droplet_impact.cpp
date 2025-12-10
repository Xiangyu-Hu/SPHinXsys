#include "2d_droplet_impact.h"
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    std::cout << "Reynolds number: " << Re_water << std::endl;
    std::cout << "Weber number: " << We_water << std::endl;
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape();
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

    FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
    air_block.defineBodyLevelSetShape();
    air_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_a, c_f), mu_a);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? air_block.generateParticles<BaseParticles, Reload>(air_block.getName())
        : air_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation water_inner(water_block);
        InnerRelation air_inner(air_block);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> water_random_particles(water_block);
        SimpleDynamics<RandomizeParticlePosition> air_random_particles(air_block);
        RelaxationStepLevelSetCorrectionInner water_relaxation_step(water_inner);
        RelaxationStepLevelSetCorrectionInner air_relaxation_step(air_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_state(sph_system);
        ReloadParticleIO write_particle_reload_files({&water_block, &air_block});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        water_random_particles.exec(0.25);
        air_random_particles.exec(0.25);
        write_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            water_relaxation_step.exec();
            air_relaxation_step.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_state.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_inner(water_block);
    ContactRelation water_air_contact(water_block, {&air_block});
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    InnerRelation air_inner(air_block);
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_wall_contact(air_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_air_complex(water_inner, {&water_air_contact, &water_wall_contact});
    ComplexRelation air_water_complex(air_inner, {&air_water_contact, &air_wall_contact});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<InitialCondition> droplet_initial_velocity(water_block);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> water_pressure_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> water_density_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> air_pressure_relaxation(air_inner, air_water_contact, air_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> air_density_relaxation(air_inner, air_water_contact, air_wall_contact);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_air_density_by_summation(air_inner, air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_water_density_by_summation(water_inner, water_air_contact, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(air_inner, air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        water_transport_correction(water_inner, water_air_contact, water_wall_contact);

    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> water_viscous_force(water_inner, water_air_contact, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> air_viscous_force(air_inner, air_water_contact, air_wall_contact);

    InteractionDynamics<fluid_dynamics::SurfaceTensionStress> water_surface_tension_stress(water_air_contact, surface_tension_coeff);
    InteractionDynamics<fluid_dynamics::SurfaceTensionStress> air_surface_tension_stress(air_water_contact, surface_tension_coeff);
    InteractionWithUpdate<fluid_dynamics::SurfaceStressForceComplex> water_surface_tension_force(water_inner, water_air_contact);
    InteractionWithUpdate<fluid_dynamics::SurfaceStressForceComplex> air_surface_tension_force(air_inner, air_water_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_water_advection_time_step_size(water_block, U_ref, 0.1);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_air_advection_time_step_size(air_block, U_ref, 0.1);
    ReduceDynamics<fluid_dynamics::SurfaceTensionTimeStep> get_water_time_step_size(water_block, 0.6);
    ReduceDynamics<fluid_dynamics::SurfaceTensionTimeStep> get_air_time_step_size(air_block, 0.6);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Matd>(water_block, "SurfaceTensionStress");
    body_states_recording.addToWrite<Vecd>(air_block, "ForcePrior");
    body_states_recording.addToWrite<Real>(air_block, "Density");
    body_states_recording.addToWrite<Real>(air_block, "Pressure");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_water_kinetic_energy(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    droplet_initial_velocity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.2;
    Real output_interval = end_time / 50; /**< Time stamps for output of body states. */
    Real dt = 0.0;                        /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_water_kinetic_energy.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            time_instance = TickCount::now();

            Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real Dt = SMIN(Dt_f, Dt_a);

            update_air_density_by_summation.exec();
            update_water_density_by_summation.exec();
            air_transport_correction.exec();
            water_transport_correction.exec();

            air_viscous_force.exec();
            water_viscous_force.exec();

            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                water_surface_tension_stress.exec();
                air_surface_tension_stress.exec();
                water_surface_tension_force.exec();
                air_surface_tension_force.exec();

                Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);

                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();

            water_block.updateCellLinkedList();
            water_air_complex.updateConfiguration();
            water_wall_contact.updateConfiguration();

            air_block.updateCellLinkedList();
            air_water_complex.updateConfiguration();
            air_wall_contact.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
        }
        write_water_kinetic_energy.writeToFile(number_of_iterations);
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_water_kinetic_energy.generateDataBase(1.0e-2);
    }
    else
    {
        write_water_kinetic_energy.testResult();
    }

    return 0;
}