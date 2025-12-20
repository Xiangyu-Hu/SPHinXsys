#include "35.h"
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
    sph_system.setReloadParticles(true);

    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));

    std::cout << "mu_f=" << mu_f << std::endl;
    std::cout << "water_block.defineBodyLevelSetShape starts" << std::endl;
    water_block.defineBodyLevelSetShape()->correctLevelSetSign();
    std::cout << "water_block.defineBodyLevelSetShape ends" << std::endl;

    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<BaseParticles, Reload>(inlet_particle_buffer, water_block.getName())
        : water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundaryFromSTL>("WallFromSTL"));

    std::cout << "wall_boundary.defineBodyLevelSetShape starts" << std::endl;
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign();
    std::cout << "wall_boundary.defineBodyLevelSetShape ends" << std::endl;

    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer_center_point(sph_system, "ObserverCenterPoint");
    observer_center_point.generateParticles<ObserverParticles>(observer_location_center_point);

    ObserverBody observer_body(sph_system, makeShared<WaterBlock>("ObserverBody")); //% Average
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? observer_body.generateParticles<BaseParticles, Reload>(water_block.getName())
        : observer_body.generateParticles<BaseParticles, Lattice>();

    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation observer_centerpoint_contact(observer_center_point, {&water_block});

    ContactRelation fluid_observer_contact2(observer_body, {&water_block}); //% Average

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
        using namespace relax_dynamics;
        /** body topology only for particle relaxation */
        InnerRelation wall_boundary_inner(wall_boundary);

        AlignedBoxByCell inlet_1_detection_box(wall_boundary,
                                                   AlignedBox(yAxis, Transform(Rotation3d(inlet_1_rotation), Vec3d(inlet_1_sub_buffer_translation)), inlet_buffer_halfsize));
        AlignedBoxByCell outlet_detection_box(wall_boundary,
                                                  AlignedBox(yAxis, Transform(Rotation3d(outlet_rotation), Vec3d(outlet_sub_buffer_translation)), outlet_buffer_halfsize));
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_boundary);
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles_water(water_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(wall_boundary);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp_water(water_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        ReloadParticleIO write_particle_reload_files_water(water_block);
        /** A  Physics relaxation step. */
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner_water(water_block_inner);

        ParticleSorting particle_sorting_wall(wall_boundary);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        random_inserted_body_particles_water.exec(0.25);

        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_inner_water.SurfaceBounding().exec();

        write_inserted_body_to_vtp.writeToFile(0);
        write_inserted_body_to_vtp_water.writeToFile(0);

        int ite_p = 0;
        int ite_max_step = 1000;
        while (ite_p < ite_max_step)
        {
            relaxation_step_inner.exec();
            relaxation_step_inner_water.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
                write_inserted_body_to_vtp_water.writeToFile(ite_p);
            }
        }

        //inlet_1_particles_deletion.exec();
        particle_sorting_wall.exec();
        wall_boundary.updateCellLinkedList();

        //outlet_particles_deletion.exec();
        particle_sorting_wall.exec();
        wall_boundary.updateCellLinkedList();

        write_inserted_body_to_vtp.writeToFile((ite_max_step + 200));
        write_inserted_body_to_vtp_water.writeToFile((ite_max_step + 200));

        std::cout << "The physics relaxation process of the wall_boundary finish !" << std::endl;
        std::cout << "The physics relaxation process of the water_block finish !" << std::endl;

        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        write_particle_reload_files_water.writeToFile(0);
        return 0;
    }

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<fluid_dynamics::DistanceFromWall> distance_to_wall(water_wall_contact);
    /** For pressure outlet . */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    /** Turbulent standard wall function needs normal vectors of wall. */
    //NearShapeSurface near_surface(water_block, makeShared<WallBoundary>("Wall"));

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_wall_contact);

    /** Pressure relaxation algorithm with Riemann solver for viscous flows. */
    //Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionForOpenBoundaryFlowWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);

    /** Density relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall<DissipativeRiemannSolver>> density_relaxation(water_block_inner, water_wall_contact);

    /** Choose one, ordinary or turbulent. Computing viscous force, */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);

    /** Impose transport velocity, with or without limiter . */
    //InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //      InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TVC_Limited_RKGC_OBC<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TVC_RKGC_OBC<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionCorrectedForOpenBoundaryFlowComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    
    /** Evaluation of density by summation approach. */
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);

    /** Initialize particle acceleration. */
    // StartupAcceleration time_dependent_acceleration(Vecd(0.0, 0.0, U_f), t_ref);
    // SimpleDynamics<GravityForce<StartupAcceleration>> apply_gravity_force(water_block, time_dependent_acceleration);

    //SimpleDynamics<InitialiseColorIndicator> initialise_color_indicator(water_block);
    //getInitialBoundingBox();
    //SimpleDynamics<InitialiseColorIndicator2> initialise_color_indicator2(water_block, box_initial_bounding);
    //----------------------------------------------------------------------
    // Inlet buffers
    //----------------------------------------------------------------------
    AlignedBox buffer_1_shape(yAxis, Transform(Rotation3d(inlet_1_rotation), Vec3d(inlet_1_buffer_translation)), inlet_buffer_halfsize);
    AlignedBoxByCell buffer_1(water_block, buffer_1_shape);
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> bidirection_buffer_1(buffer_1, inlet_particle_buffer);
    //SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> buffer_1_inflow_pressure_condition(buffer_1);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<LeftInflowPressure>> buffer_1_inflow_pressure_condition(buffer_1);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> buffer_1_inflow_velocity_condition(buffer_1);
    //----------------------------------------------------------------------
    // Outlet buffer
    //----------------------------------------------------------------------
    AlignedBox outlet_buffer_shape(yAxis, Transform(Rotation3d(outlet_rotation_reverse), Vecd(outlet_buffer_translation)), outlet_buffer_halfsize);
    AlignedBoxByCell outlet_buffer(water_block, outlet_buffer_shape);
    fluid_dynamics::BidirectionalBuffer<RightOutflowPressure> outlet_bidirection_buffer(outlet_buffer, inlet_particle_buffer);

    //SimpleDynamics<fluid_dynamics::PressureCondition<RightOutflowPressure>> outflow_pressure_condition(outlet_buffer);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<RightOutflowPressure>> outflow_pressure_condition(outlet_buffer);
    //----------------------------------------------------------------------

    //SimpleDynamics<ClearBufferParticleIndicator> clear_buffer_particle_indicator(water_block, zAxis, H_inlet, H_inlet + BW); //% This is case-dependent

    //InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density_pressure(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density_freestream(water_block_inner, water_wall_contact);

    /** Choose one, ordinary or turbulent. Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);

    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);

    ObservingAQuantity<Real> observing_pressure(fluid_observer_contact2, "Pressure");          //% Average pressure
    SimpleDynamics<ParticleSnapshotAverage<Real>> average_pressure(observer_body, "Pressure"); //% Average pressure
    //ObservingAQuantity<int> observing_buffer_particle_indicator(fluid_observer_contact2, "BufferIndicator");          //% Average

    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    ParticleSorting particle_sorting_wall(wall_boundary);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure"); // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator"); // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Density");  // output for debug
    //body_states_recording.addToWrite<Vecd>(water_block, "ZeroGradientResidue"); // output for debug
    // ObservedQuantityRecording<Vecd> write_recorded_water_centerline_velocity("Velocity", fluid_observer_centerline_contact);
    // ObservedQuantityRecording<Real> write_recorded_water_centerline_pressure("Pressure", fluid_observer_centerline_contact);
    // ObservedQuantityRecording<Vecd> write_recorded_water_velocity_cross_section("Velocity", fluid_observer_cross_section_contact);
    body_states_recording.addToWrite<int>(water_block, "BufferIndicator");
    //body_states_recording.addToWrite<Real>(water_block, "VolumetricMeasure");
    body_states_recording.addToWrite<Matd>(water_block, "LinearGradientCorrectionMatrix");

    WriteToVtpIfVelocityOutOfBound abnormal_velocity_recording(sph_system, 1.0e6 * U_max);

    //% Temporary Treat Note that this should be in front of TAG particle include inlet outlet and buffer
    //SimpleDynamics<DisposerForInitialParticleDeletion> delete_initial_particle(water_block);
    //delete_initial_particle.exec();
    //particle_sorting.exec();
    //water_block.updateCellLinkedList();

    //initialise_color_indicator.exec();
    //initialise_color_indicator2.exec();

    BodyStatesRecordingToVtp write_observation_states(observer_body);     //% Average
    write_observation_states.addToWrite<Real>(observer_body, "Pressure"); //% Average pressure
    //write_observation_states.addToWrite<int>(observer_body, "BufferIndicator");  //% Average
    // body_states_recording.addToWrite<int>(observer_body, "BufferIndicator"); //% Average

    /**
     * @brief Setup geometry and initial conditions.
     */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");

    /** Tag inlet/outlet truncated particles */
    inlet_outlet_surface_particle_indicator.exec();
    /** Tag in/outlet buffer particles */
    bidirection_buffer_1.tag_buffer_particles.exec();

    outlet_bidirection_buffer.tag_buffer_particles.exec();

    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 6000.0;                     /**< End time. */
    Real cutoff_ratio = 0.92;                   //** cutoff_time should be a integral and the same as the PY script */
    Real cutoff_time = cutoff_ratio * end_time; //** cutoff_time should be a integral and the same as the PY script */
    Real num_output_files = 600000.0;
    Real Output_Time = end_time / num_output_files; /**< Time stamps for output of body states. */
    Real index_check_file_fully_developed = num_output_files * cutoff_ratio;
    Real dt = 0.0; /**< Default acoustic time step sizes. */

    Real time_output_average_data = 400.0; //% Average

    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    particle_sorting_wall.exec();
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
     //std::cout << "Simulation starts?" << std::endl;
     //std::cin.get();
    int num_output_file = 0;
    std::ofstream logfile("output/output.log");
    std::ofstream mixing_file("output/mixing_rate.dat");
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            //apply_gravity_force.exec();

            Real Dt = get_fluid_advection_time_step_size.exec();
            //Real Dt = get_turbulent_fluid_advection_time_step_size.exec();

            //update_density_by_summation.exec();
            //update_fluid_density_pressure.exec();
            update_fluid_density_freestream.exec();

            corrected_configuration_fluid.exec();

            viscous_force.exec();
            //turbulent_viscous_force.exec();

            transport_velocity_correction.exec();

            kernel_summation.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            int inner_itr = 0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                pressure_relaxation.exec(dt);

                buffer_1_inflow_pressure_condition.exec(dt);

                outflow_pressure_condition.exec(dt);

                buffer_1_inflow_velocity_condition.exec();

                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_itr++;
                //std::cout << "num_output_file=" << num_output_file << std::endl;
                //if (physical_time >9.3)
                //{
                //body_states_recording.writeToFile();
                //}
                abnormal_velocity_recording.writeToFile();
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                logfile << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                        << physical_time
                        << "	Dt = " << Dt << "	dt = " << dt << std::endl;
            }
            number_of_iterations++;

            // ** First do injection for all buffers *
            bidirection_buffer_1.injection.exec();
            outlet_bidirection_buffer.injection.exec();
            // ** Then do deletion for all buffers *
            bidirection_buffer_1.deletion.exec();
            outlet_bidirection_buffer.deletion.exec();

            /** Update cell linked list and configuration. */
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                //particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            //fluid_observer_centerline_contact.updateConfiguration();
            //fluid_observer_cross_section_contact.updateConfiguration();

            /** Tag truncated inlet/outlet particles*/
            inlet_outlet_surface_particle_indicator.exec();

            //% Clear before tagging.
            //clear_buffer_particle_indicator.exec();

            /** Tag in/outlet buffer particles that suffer pressure condition*/
            bidirection_buffer_1.tag_buffer_particles.exec();
            outlet_bidirection_buffer.tag_buffer_particles.exec();

            if (physical_time > time_output_average_data * 100.0)
            {
                fluid_observer_contact2.updateConfiguration(); //% Average
                //% Average pressure
                observing_pressure.exec();
                average_pressure.exec();
                // observing_buffer_particle_indicator.exec();
            }

            // if (physical_time > cutoff_time)
            // {
            //     write_recorded_water_centerline_velocity.writeToFile(number_of_iterations);
            //     write_recorded_water_centerline_pressure.writeToFile(number_of_iterations);
            //     write_recorded_water_velocity_cross_section.writeToFile(number_of_iterations);
            // }
            //if (physical_time > end_time * 0.5)
            //body_states_recording.writeToFile();
        }
        //TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        observer_centerpoint_contact.updateConfiguration();
        num_output_file++;
        //if (num_output_file == 100)
        //    system("pause");
        //TickCount t3 = TickCount::now();
        
        if (physical_time > time_output_average_data)
        {
            write_observation_states.writeToFile(); //% Average
        }
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << "Cutoff_time: " << cutoff_time
              << " seconds." << std::endl;
    std::cout << "For checking fully-developed or not, index of the cutoff output file =  " << index_check_file_fully_developed << std::endl;
    logfile << "Total wall time for computation: " << tt.seconds()
            << " seconds." << std::endl;
    logfile.close();
    mixing_file.close();
    return 0;
}
