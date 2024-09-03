#include "2d_turbulent_channel_VI_offset_model.h"
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

    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    IOEnvironment io_environment(sph_system);
    /**
     * @brief Material property, particles and body creation of fluid.
     */

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape();
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticlesWithReserve<BaseParticles, Reload>(inlet_particle_buffer, water_block.getName())
        : water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>( wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer_center_point(sph_system, "ObserverCenterPoint");
    observer_center_point.generateParticles<ObserverParticles>(observer_location_center_point);

   for(int i = 0; i < num_observer_points; ++i)
   {
       Vecd pos_observer_i= pos_observe_start + i * observe_spacing * unit_direction_observe;
       if( i==0 )
       {
           pos_observer_i -= observer_offset_distance * unit_direction_observe ;
       }
       if( i == num_observer_points - 1 )
       {
           pos_observer_i += observer_offset_distance * unit_direction_observe ;
       }
       observation_location.push_back(pos_observer_i);
   }
    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_location);


    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    ContactRelation observer_centerpoint_contact(observer_center_point, {&water_block});
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
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_boundary);
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles_water(water_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp( wall_boundary );
        BodyStatesRecordingToVtp write_inserted_body_to_vtp_water(water_block );
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files( wall_boundary);
        ReloadParticleIO write_particle_reload_files_water(water_block);
        /** A  Physics relaxation step. */
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner_water(water_block_inner);
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
        while (ite_p < 1000)
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
    //InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(ConstructorArgs(water_block_inner, 0.9), water_wall_contact);
    
    //InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_fluid(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::TurbulentLinearGradientCorrectionMatrixInner> corrected_configuration_fluid_only_inner(water_block_inner);

    /** Pressure relaxation algorithm with Riemann solver for viscous flows. */
    //Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionForOpenBoundaryFlowWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    
    /** Density relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);    
    //Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_wall_contact);    
    //Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall<DissipativeRiemannSolver>> density_relaxation(water_block_inner, water_wall_contact);    
    
    
    /** Turbulent.Note: When use wall function, K Epsilon calculation only consider inner */
    InteractionWithUpdate<fluid_dynamics::JudgeIsNearWall> update_near_wall_status(water_block_inner, water_wall_contact);

    InteractionWithUpdate<fluid_dynamics::GetVelocityGradientInner> get_velocity_gradient(water_block_inner);
    //InteractionWithUpdate<fluid_dynamics::GetVelocityGradientComplex> get_velocity_gradient(water_block_inner, water_wall_contact);

    //InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<LinearGradientCorrection>> vel_grad_calculation(water_block_inner, water_wall_contact);
    //SimpleDynamics<fluid_dynamics::TransferVelocityGradient> transfer_velocity_gradient(water_block);


    /** Turbulent.Note: Temporarily transfer parameters at this place. The 3rd parameter refers to extra dissipation for viscous */
    
    InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelInner> k_equation_relaxation(water_block_inner, initial_turbu_values);
    //InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelInner> k_equation_relaxation(water_block_inner, initial_turbu_values, 1);

    InteractionWithUpdate<fluid_dynamics::E_TurtbulentModelInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::TKEnergyForceComplex> turbulent_kinetic_energy_force(water_block_inner, water_wall_contact);
    InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection> standard_wall_function_correction(water_block_inner, water_wall_contact, y_p_constant);

    SimpleDynamics<fluid_dynamics::ConstrainNormalVelocityInRegionP> constrain_normal_velocity_in_P_region(water_block);
    //SimpleDynamics<fluid_dynamics::ConstrainVelocityAt_Y_Direction> constrain_Y_velocity(water_block, DL);
    //SimpleDynamics<fluid_dynamics::UpdateTurbulentPlugFlowIndicator> update_turbu_plug_flow_indicator(water_block, DH);

    SimpleDynamics<fluid_dynamics::GetTimeAverageCrossSectionData> get_time_average_cross_section_data(water_block_inner, num_observer_points, monitoring_bound,offset_distance);

    /** Choose one, ordinary or turbulent. Computing viscous force, */
    InteractionWithUpdate<fluid_dynamics::TurbulentViscousForceWithWall> turbulent_viscous_force(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);
    
    /** Impose transport velocity, with or without Extra Transprot Force. */
    //InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::ExtraTransportForceComplex<BulkParticles>> impose_extra_transport_force(water_block_inner, water_wall_contact);

    /** Impose transport velocity, with or without Extra Transprot Force, and with limiter . */
    //InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::ExtraTransportForceLimitedComplex<BulkParticles>> impose_extra_transport_force(water_block_inner, water_wall_contact);

    //InteractionWithUpdate<fluid_dynamics::TVC_Limited_withLinearGradientCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TVC_NoLimiter_withLinearGradientCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TVC_ModifiedLimited_withLinearGradientCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TVC_ModifiedLimited_withoutLinearGradientCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TVC_ModifiedLimited_RKGC_OBFCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TVC_NotLimited_RKGC_OBFCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::TVC_ModifiedLimited_NoRKGC<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);

    /** A temporarily test for the limiter . */
    SimpleDynamics<fluid_dynamics::GetLimiterOfTransportVelocityCorrection> get_limiter_of_transport_velocity_correction(water_block, 1000);
    

    /** Evaluation of density by summation approach. */
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);
    
    /** Initialize particle acceleration. */
    TimeDependentAcceleration time_dependent_acceleration(Vec2d::Zero());
    SimpleDynamics<GravityForce> apply_gravity_force(water_block, time_dependent_acceleration);
    
    //----------------------------------------------------------------------
    // Left/Inlet buffer
    //----------------------------------------------------------------------
    //BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    //SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer, xAxis);
    //BodyAlignedBoxByCell inlet_velcoity_buffer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    //SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inlet_velocity_buffer_inflow_condition(inlet_velcoity_buffer, 1.0);
    
    BodyAlignedBoxByCell left_emitter(water_block, makeShared<AlignedBoxShape>(xAxis,Transform(Vec2d(left_buffer_translation)), left_buffer_halfsize));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer left_emitter_inflow_injection(left_emitter, inlet_particle_buffer);

    BodyAlignedBoxByCell left_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(Pi), Vec2d(left_buffer_translation)), left_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer);
    
    //SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);

    /** Turbulent InflowTurbulentCondition.It needs characteristic Length to calculate turbulent length  */
    SimpleDynamics<fluid_dynamics::InflowTurbulentCondition> impose_turbulent_inflow_condition(left_emitter, characteristic_length, 0.8);
    
    //----------------------------------------------------------------------
    // Right/Outlet buffer
    //----------------------------------------------------------------------
    //BodyAlignedBoxByCell disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    //SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, 0);
    BodyAlignedBoxByCell right_emitter(water_block, makeShared<AlignedBoxShape>(xAxis,Transform(Rotation2d(Pi), Vec2d(right_buffer_translation)), right_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightOutflowPressure> right_emitter_inflow_injection(right_emitter, inlet_particle_buffer);
    BodyAlignedBoxByCell right_disposer(water_block, makeShared<AlignedBoxShape>(xAxis,Transform(Vec2d(right_buffer_translation)), right_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer);
    
    //SimpleDynamics<fluid_dynamics::PressureCondition<RightOutflowPressure>> right_outflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<RightOutflowPressure>> right_outflow_pressure_condition(right_emitter);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density_pressure(water_block_inner, water_wall_contact);

    /** Choose one, ordinary or turbulent. Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    //ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity> update_eddy_viscosity(water_block);

    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");		   // output for debug
    body_states_recording.addToWrite<int>(water_block,"Indicator"); // output for debug
    body_states_recording.addToWrite<Real>(water_block,"Density"); // output for debug
    body_states_recording.addToWrite<Vecd>(water_block,"ZeroGradientResidue"); // output for debug
    ObservedQuantityRecording<Vecd> write_recorded_water_velocity("Velocity", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_k("TurbulenceKineticEnergy", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_mut("TurbulentViscosity", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_epsilon("TurbulentDissipation", fluid_observer_contact);
    body_states_recording.addToWrite<int>(water_block, "BufferParticleIndicator");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_centerpoint_quantity("TurbulentViscosity", observer_centerpoint_contact);
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
    left_emitter_inflow_injection.tag_buffer_particles.exec();
    right_emitter_inflow_injection.tag_buffer_particles.exec();

    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 200.0;   /**< End time. */
    Real Output_Time = end_time / 4.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_centerpoint_quantity.writeToFile(number_of_iterations);
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
            apply_gravity_force.exec();

            //Real Dt = get_fluid_advection_time_step_size.exec();
            Real Dt = get_turbulent_fluid_advection_time_step_size.exec();

            //inlet_outlet_surface_particle_indicator.exec();

            //update_density_by_summation.exec();
            update_fluid_density_pressure.exec();

            corrected_configuration_fluid.exec();
            corrected_configuration_fluid_only_inner.exec();

            update_eddy_viscosity.exec();

            //viscous_force.exec();
            turbulent_viscous_force.exec();

            transport_velocity_correction.exec();
            get_limiter_of_transport_velocity_correction.exec();

            //impose_extra_transport_force.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            int inner_itr = 0;
            while (relaxation_time < Dt)
            {
                //if (inner_itr == 8)
                //{
                //    system("pause");
                //}

                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                turbulent_kinetic_energy_force.exec();

                pressure_relaxation.exec(dt);

                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_outflow_pressure_condition.exec(dt);

                //constrain_normal_velocity_in_P_region.exec();
                
                //** For test *
                //constrain_Y_velocity.exec();
                //update_turbu_plug_flow_indicator.exec();

                //inlet_velocity_buffer_inflow_condition.exec();
                inflow_velocity_condition.exec();

                impose_turbulent_inflow_condition.exec();

                density_relaxation.exec(dt);

                distance_to_wall.exec();
                update_near_wall_status.exec();

                standard_wall_function_correction.exec();
                
                get_velocity_gradient.exec(dt);
                //vel_grad_calculation.exec();
                //** To test the latest vel_grad without changing the current code framework */
                //transfer_velocity_gradient.exec();
                
                k_equation_relaxation.exec(dt);
                epsilon_equation_relaxation.exec(dt);
                //k_equation_relaxation.update_prior_turbulent_value();

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                inner_itr++;
                //std::cout << "num_output_file=" << num_output_file << std::endl;
                //if (GlobalStaticVariables::physical_time_ >9.3)
                //{
                    //body_states_recording.writeToFile();
                //}
                
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerpoint_quantity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            /** inflow injection*/
            //emitter_inflow_injection.exec();
            left_emitter_inflow_injection.injection.exec();
            right_emitter_inflow_injection.injection.exec();

            //disposer_outflow_deletion.exec();
            left_disposer_outflow_deletion.exec();
            right_disposer_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();

            /** Tag truncated inlet/outlet particles*/
            inlet_outlet_surface_particle_indicator.exec();
            /** Tag in/outlet buffer particles that suffer pressure condition*/
            left_emitter_inflow_injection.tag_buffer_particles.exec();
            right_emitter_inflow_injection.tag_buffer_particles.exec();

            if (GlobalStaticVariables::physical_time_ > end_time * 0.6) 
            {
                get_time_average_cross_section_data.exec();
                get_time_average_cross_section_data.output_time_history_data(end_time * 0.75);
                write_recorded_water_velocity.writeToFile(number_of_iterations);
                write_recorded_water_k.writeToFile(number_of_iterations);
                write_recorded_water_mut.writeToFile(number_of_iterations);
                write_recorded_water_epsilon.writeToFile(number_of_iterations);
            }
            //if (GlobalStaticVariables::physical_time_ > end_time * 0.5)
                //body_states_recording.writeToFile();


        }
        //TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        observer_centerpoint_contact.updateConfiguration();
        num_output_file++;
        //if (num_output_file == 100)
        //    system("pause");
        //TickCount t3 = TickCount::now();

    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    get_time_average_cross_section_data.get_time_average_data(end_time * 0.75);
    std::cout << "The time-average data is output " << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_centerpoint_quantity.generateDataBase(1.0e-3);
    }
    else
    {
        write_centerpoint_quantity.testResult();
    }
    return 0;
}
