#include "2d_turbulent_strongly_curved_channel.h"
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
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    getPositionsOfMultipleObserveLines();
    output_observe_positions();
    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_locations);

    /** topology */
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
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(wall_boundary);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp_water(water_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        ReloadParticleIO write_particle_reload_files_water(water_block);
        /** A  Physics relaxation step. */
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner_water(water_block_inner);
        //RelaxationStepLevelSetCorrectionComplex relaxation_step_inner(wall_boundary_inner, water_block);
        //RelaxationStepLevelSetCorrectionComplex relaxation_step_inner_water(water_block_inner, water_wall_contact); //** For consistency of Fluid *

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

    InteractionWithUpdate<fluid_dynamics::TurbulentLinearGradientCorrectionMatrixInner> corrected_configuration_fluid_only_inner(water_block_inner);

    /** Pressure relaxation algorithm with Riemann solver for viscous flows. */
    //Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    //Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionForOpenBoundaryFlowWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);

    /** Density relaxation algorithm by using position verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);

    /** Turbulent.Note: When use wall function, K Epsilon calculation only consider inner */
    InteractionWithUpdate<fluid_dynamics::JudgeIsNearWall> update_near_wall_status(water_block_inner, water_wall_contact);

    InteractionWithUpdate<fluid_dynamics::GetVelocityGradientInner> get_velocity_gradient(water_block_inner, weight_vel_grad_sub_nearwall);
    //InteractionWithUpdate<fluid_dynamics::GetVelocityGradientComplex> get_velocity_gradient(water_block_inner, water_wall_contact);

    /** Turbulent.Note: Temporarily transfer parameters at this place. The 3rd parameter refers to extra dissipation for viscous */
    InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelInner> k_equation_relaxation(water_block_inner, initial_turbu_values, is_AMRD);
    InteractionWithUpdate<fluid_dynamics::E_TurtbulentModelInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::TKEnergyForceComplex> turbulent_kinetic_energy_force(water_block_inner, water_wall_contact);
    InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection> standard_wall_function_correction(water_block_inner, water_wall_contact, y_p_constant);

    /** Turbulent.Extra boundary condition */
    //SimpleDynamics<fluid_dynamics::ConstrainNormalVelocityInRegionP> constrain_normal_velocity_in_P_region(water_block);

    /** Choose one, ordinary or turbulent. Computing viscous force, */
    InteractionWithUpdate<fluid_dynamics::TurbulentViscousForceWithWall> turbulent_viscous_force(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);

    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TVC_ModifiedLimited_RKGC_OBFCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);

    /** A temporarily test for the limiter . */
    SimpleDynamics<fluid_dynamics::GetLimiterOfTransportVelocityCorrection> get_limiter_of_transport_velocity_correction(water_block, 1000);

    /** Initialize particle acceleration. */
    TimeDependentAcceleration time_dependent_acceleration(Vec2d::Zero());
    SimpleDynamics<GravityForce> apply_gravity_force(water_block, time_dependent_acceleration);

    //----------------------------------------------------------------------
    // Left/Inlet buffer
    //----------------------------------------------------------------------
    BodyAlignedBoxByCell left_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(left_buffer_translation)), left_buffer_halfsize));
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer left_emitter_inflow_injection(left_emitter, inlet_particle_buffer);

    BodyAlignedBoxByCell left_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(Pi), Vec2d(left_buffer_translation)), left_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer);

    //SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);

    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);

    /** Turbulent InflowTurbulentCondition.It needs characteristic Length to calculate turbulent length  */
    SimpleDynamics<fluid_dynamics::InflowTurbulentCondition> impose_turbulent_inflow_condition(left_emitter, characteristic_length, relaxation_rate_turbulent_inlet, type_turbulent_inlet);

    //----------------------------------------------------------------------
    // Right/Outlet buffer
    //----------------------------------------------------------------------
    //BodyAlignedBoxByCell disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    //SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, 0);
    BodyAlignedBoxByCell right_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(outlet_emitter_rotation_angel), Vec2d(right_buffer_translation)), right_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightOutflowPressure> right_emitter_inflow_injection(right_emitter, inlet_particle_buffer);
    BodyAlignedBoxByCell right_disposer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(outlet_disposer_rotation_angel), Vec2d(right_buffer_translation)), right_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer);

    //SimpleDynamics<fluid_dynamics::PressureCondition<RightOutflowPressure>> right_outflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<RightOutflowPressure>> right_outflow_pressure_condition(right_emitter);

    /** Temporary treatment for Pressure outlet module  */
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density_pressure(water_block_inner, water_wall_contact);
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);

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
    body_states_recording.addToWrite<Real>(water_block, "Pressure");            // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");            // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Density");             // output for debug
    body_states_recording.addToWrite<Vecd>(water_block, "ZeroGradientResidue"); // output for debug
    ObservedQuantityRecording<Vecd> write_recorded_water_velocity("Velocity", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_k("TurbulenceKineticEnergy", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_mut("TurbulentViscosity", fluid_observer_contact);
    ObservedQuantityRecording<Real> write_recorded_water_epsilon("TurbulentDissipation", fluid_observer_contact);
    body_states_recording.addToWrite<int>(water_block, "BufferParticleIndicator");

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

    /** Output the start states of bodies. */
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    Real dt = 0.0; /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    int num_output_file = 0;
    //Real start_time_turbulence = 70.0;
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

                inflow_velocity_condition.exec();

                impose_turbulent_inflow_condition.exec();

                density_relaxation.exec(dt);

                distance_to_wall.exec();
                update_near_wall_status.exec();

                if (GlobalStaticVariables::physical_time_ > turbulent_module_activate_time) //** A temporary treatment *
                {
                    standard_wall_function_correction.exec();
                    get_velocity_gradient.exec(dt);
                    k_equation_relaxation.exec(dt);
                    epsilon_equation_relaxation.exec(dt);
                }
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
            }
            number_of_iterations++;

            /** inflow injection*/
            left_emitter_inflow_injection.injection.exec();
            right_emitter_inflow_injection.injection.exec();

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

            if (GlobalStaticVariables::physical_time_ > cutoff_time)
            {
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
    return 0;
}
