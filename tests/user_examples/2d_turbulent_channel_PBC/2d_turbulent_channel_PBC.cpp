#include "2d_turbulent_channel_PBC.h"
using namespace SPH;

int main(int ac, char *av[])
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem sph_system(system_domain_bounds, resolution_ref);

    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);

    sph_system.setGenerateRegressionData(false);

    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    IOEnvironment io_environment(sph_system);
    /**
     * @brief Material property, particles and body creation of fluid.
     */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineComponentLevelSetShape("OuterBoundary");
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<ParticleGeneratorReload>(wall_boundary.getName())
        : wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd(DL / 2.0, DH / 2.0)};
    fluid_observer.generateParticles<ParticleGeneratorObserver>(observation_location);

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        //InnerRelation water_block_inner_relax(water_block);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(water_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp({&water_block});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(water_block);
        /** A  Physics relaxation step. */
        //relax_dynamics::RelaxationStepInner relaxation_step_inner(water_block_inner_relax);
        relax_dynamics::RelaxationStepLevelSetCorrectionComplex relaxation_step_complex(
            ConstructorArgs(water_block_inner, "OuterBoundary"), water_wall_contact);

        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_complex.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_complex.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of the water_block finish !" << std::endl;

        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }

    /**
     * @brief 	Define all numerical methods which are used in this case.
    /** Pressure relaxation algorithm with Riemann solver for viscous flows. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    /** Density relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    /** Turbulent standard wall function needs normal vectors of wall. */
    NearShapeSurface near_surface(water_block, makeShared<WallBoundary>("Wall"));
    near_surface.level_set_shape_.writeLevelSet(sph_system);

    /** Turbulent.Note: When use wall function, K Epsilon calculation only consider inner */
    InteractionWithUpdate<fluid_dynamics::JudgeIsNearWall> update_near_wall_status(water_block_inner, water_wall_contact, near_surface);
    InteractionDynamics<fluid_dynamics::GetVelocityGradientInner> get_velocity_gradient(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::K_TurbulentModelInner> k_equation_relaxation(water_block_inner, initial_turbu_values);
    InteractionWithUpdate<fluid_dynamics::E_TurbulentModelInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::TKEnergyForceComplex> turbulent_kinetic_energy_force(water_block_inner, water_wall_contact);
    SimpleDynamics<fluid_dynamics::StandardWallFunctionCorrection> standard_wall_function_correction(water_block, offset_dist_ref);

    SimpleDynamics<fluid_dynamics::GetTimeAverageCrossSectionData, SequencedPolicy> get_time_average_cross_section_data(water_block_inner, num_observer_points, monitoring_bound);

    /** Choose one, ordinary or turbulent. Computing viscous force, */
    InteractionWithUpdate<fluid_dynamics::TurbulentViscousForceWithWall> turbulent_viscous_force(water_block_inner, water_wall_contact);
    //InteractionDynamics<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);

    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_wall_contact);
    water_block.addBodyStateForRecording<Real>("Pressure"); // output for debug
    water_block.addBodyStateForRecording<int>("Indicator"); // output for debug
    water_block.addBodyStateForRecording<Real>("Density");  // output for debug

    /** Initialize particle force. */
    Gravity gravity(Vecd(gravity_g, 0.0));
    SimpleDynamics<GravityForce> constant_gravity(water_block, gravity);

    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, water_block.getBodyShapeBounds(), xAxis);

    /** Choose one, ordinary or turbulent. Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    //ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);

    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity> update_eddy_viscosity(water_block);

    //SimpleDynamics<fluid_dynamics::ClearYPositionForTest> clear_y_displacement(water_block);

    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);

    //RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>>
    //write_water_mechanical_energy(water_block);

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 600.0;              /**< End time. */
    Real Output_Time = end_time / 40.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                      /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    //write_water_mechanical_energy.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    int num_output_file = 0;
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            //Real Dt = get_fluid_advection_time_step_size.exec();
            Real Dt = get_turbulent_fluid_advection_time_step_size.exec();

            update_density_by_summation.exec();

            update_eddy_viscosity.exec();

            //viscous_force.exec();
            turbulent_viscous_force.exec();

            transport_velocity_correction.exec();
            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            int inner_itr = 0;
            while (relaxation_time < Dt)
            {
                //std::cout << "inner_itr=" << inner_itr << std::endl;

                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                turbulent_kinetic_energy_force.exec();

                pressure_relaxation.exec(dt);

                density_relaxation.exec(dt);

                update_near_wall_status.exec();
                get_velocity_gradient.exec(dt);
                standard_wall_function_correction.exec();
                k_equation_relaxation.exec(dt);
                epsilon_equation_relaxation.exec(dt);
                k_equation_relaxation.update_prior_turbulent_value();

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                inner_itr++;

                //clear_y_displacement.exec();
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                //if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                //{
                // write_water_mechanical_energy.writeToFile(number_of_iterations);
                //}
            }
            number_of_iterations++;
            /** Update cell linked list and configuration. */
            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            periodic_condition.update_cell_linked_list_.exec();
            water_block_complex.updateConfiguration();

            get_time_average_cross_section_data.exec();
            get_time_average_cross_section_data.output_time_history_data(end_time * 0.75);
            fluid_observer_contact.updateConfiguration();
        }
        TickCount t2 = TickCount::now();

        body_states_recording.writeToFile();
        num_output_file++;
        if (num_output_file == 20)
            //system("pause");

            TickCount t3 = TickCount::now();
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    get_time_average_cross_section_data.get_time_average_data(end_time * 0.75);
    std::cout << "The time-average data is output " << std::endl;

    //if (sph_system.GenerateRegressionData())
    //{
    //    write_water_mechanical_energy.generateDataBase(0.1);
    //}
    //else if (sph_system.RestartStep() == 0)
    //{
    //    write_water_mechanical_energy.testResult();
    //}

    return 0;
}
