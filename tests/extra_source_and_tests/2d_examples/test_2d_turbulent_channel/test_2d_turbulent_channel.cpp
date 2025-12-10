#include "test_2d_turbulent_channel.h"
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

    sph_system.handleCommandlineOptions(ac, av);
    /**
     * @brief Material property, particles and body creation of fluid.
     */

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape();
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
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

    ObserverBody observer_center_point(sph_system, "ObserverCenterPoint");
    observer_center_point.generateParticles<ObserverParticles>(observer_location_center_point);

    /** topology */
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
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
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(wall_boundary);
        BodyStatesRecordingToVtp write_inserted_body_to_vtp_water(water_block);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(wall_boundary);
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
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TurbulentLinearGradientCorrectionMatrixInner> corrected_configuration_fluid_only_inner(water_block_inner);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionForOpenBoundaryFlowWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);

    InteractionWithUpdate<fluid_dynamics::JudgeIsNearWall> update_near_wall_status(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::GetVelocityGradientInner> get_velocity_gradient(water_block_inner, weight_vel_grad_sub_nearwall);
    InteractionWithUpdate<fluid_dynamics::K_TurbulentModelInner> k_equation_relaxation(water_block_inner, initial_turbu_values, is_AMRD, is_source_term_linearisation);
    InteractionWithUpdate<fluid_dynamics::E_TurbulentModelInner> epsilon_equation_relaxation(water_block_inner, is_source_term_linearisation);
    InteractionDynamics<fluid_dynamics::TKEnergyForceComplex> turbulent_kinetic_energy_force(water_block_inner, water_wall_contact);
    InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection> standard_wall_function_correction(water_block_inner, water_wall_contact, y_p_constant);
    SimpleDynamics<fluid_dynamics::ConstrainNormalVelocityInRegionP> constrain_normal_velocity_in_P_region(water_block);
    InteractionWithUpdate<fluid_dynamics::TurbulentViscousForceWithWall> turbulent_viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TVC_ModifiedLimited_RKGC_OBFCorrection<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);

    StartupAcceleration time_dependent_acceleration(Vec2d(U_f, 0.0), 2.0);
    SimpleDynamics<GravityForce<StartupAcceleration>> apply_gravity_force(water_block, time_dependent_acceleration);

    //----------------------------------------------------------------------
    // Left/Inlet buffer
    //----------------------------------------------------------------------
    AlignedBox left_emitter_shape(xAxis, Transform(Vec2d(left_buffer_translation)), left_buffer_halfsize);
    AlignedBoxByCell left_emitter(water_block, left_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, inlet_particle_buffer);

    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::InflowTurbulentCondition> impose_turbulent_inflow_condition(left_emitter, characteristic_length, relaxation_rate_turbulent_inlet, type_turbulent_inlet);

    //----------------------------------------------------------------------
    // Right/Outlet buffer
    //----------------------------------------------------------------------
    AlignedBox right_emitter_shape(xAxis, Transform(Rotation2d(Pi), Vec2d(right_buffer_translation)), right_buffer_halfsize);
    AlignedBoxByCell right_emitter(water_block, right_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<RightOutflowPressure> right_bidirection_buffer(right_emitter, inlet_particle_buffer);
    SimpleDynamics<fluid_dynamics::PressureConditionCorrection<RightOutflowPressure>> right_outflow_pressure_condition(right_emitter);
    //----------------------------------------------------------------------

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density_pressure(water_block_inner, water_wall_contact);
    ReduceDynamics<fluid_dynamics::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity> update_eddy_viscosity(water_block);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	File output and regression check.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "TurbulenceKineticEnergy");
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

    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();

    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 200.0;             /**< End time. */
    Real Output_Time = end_time / 4.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;
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
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            apply_gravity_force.exec();

            Real Dt = get_turbulent_fluid_advection_time_step_size.exec();

            update_fluid_density_pressure.exec();

            corrected_configuration_fluid.exec();
            corrected_configuration_fluid_only_inner.exec();

            update_eddy_viscosity.exec();

            turbulent_viscous_force.exec();

            transport_velocity_correction.exec();

            Real relaxation_time = 0.0;
            int inner_itr = 0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                turbulent_kinetic_energy_force.exec();

                pressure_relaxation.exec(dt);

                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_outflow_pressure_condition.exec(dt);

                if (is_constrain_normal_velocity_in_P_region)
                    constrain_normal_velocity_in_P_region.exec();

                inflow_velocity_condition.exec();

                impose_turbulent_inflow_condition.exec();

                density_relaxation.exec(dt);

                distance_to_wall.exec();
                update_near_wall_status.exec();

                standard_wall_function_correction.exec();

                get_velocity_gradient.exec(dt);

                k_equation_relaxation.exec(dt);
                epsilon_equation_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_itr++;
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_centerpoint_quantity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();

            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();

            inlet_outlet_surface_particle_indicator.exec();

            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }

        body_states_recording.writeToFile();
        observer_centerpoint_contact.updateConfiguration();
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

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
