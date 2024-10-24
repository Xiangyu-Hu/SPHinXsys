/**
 * @file 	2d_turbulent_flow_around_cylinder.cpp
 * @brief 	2d_turbulent_flow_around_cylinder flow with K-Epsilon two equations RANS model.
 * @details This is the one of the basic test cases.
 * @author 	Xiangyu Hu
 */
#include "2d_turbulent_flow_around_cylinder.h"

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);

    /** Tag for run particle relaxation for the initial body fitted distribution. */
    system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    system.setReloadParticles(true);
    /** handle command line arguments. */
    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    //water_block.defineAdaptationRatios(1.4);

    SolidBody cylinder(system, makeShared<Cylinder>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 2.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? cylinder.generateParticles<ParticleGeneratorReload>(io_environment, cylinder.getName())
        : cylinder.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(system, "FluidObserver");
    fluid_observer.generateParticles<ParticleGeneratorObserver>(observation_locations);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which are only use for now for updating configurations.
    //--------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(io_environment, {&cylinder});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, {&cylinder});
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------

    /** Initialize particle acceleration. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d::Zero()));
    /** Emitter buffer inflow condition. */
    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
    BodyAlignedBoxByCell emitter_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition(emitter_buffer);
    /** Emitter buffer outflow condition. */
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, xAxis);
    /** time-space method to detect surface particles. */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density(water_block_inner, water_block_contact);
    /** Output a method-specific particle data for debug */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("Indicator");

    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** modify the velocity of boundary particles with free-stream velocity. */
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint(water_block);

    //Attention! the original one does not use Riemann solver for pressure
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);

    /** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
    pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);

    //Attention! the original one does use Riemann solver for density
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);

    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    /** compute the vorticity. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);

    /** Turbulent.Note: When use wall function, K Epsilon calculation only consider inner */
    InteractionWithUpdate<fluid_dynamics::K_TurbulentModelInner> k_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::GetVelocityGradientInner> get_velocity_gradient(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::E_TurbulentModelInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::TKEnergyAccComplex> turbulent_kinetic_energy_acceleration(water_block_inner, water_block_contact);

    /** Turbulent advection time step. */
    ReduceDynamics<fluid_dynamics::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    //ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);

    /** Turbulent standard wall function needs normal vectors of wall. */
    NearShapeSurface near_surface(water_block, makeShared<Cylinder>("Cylinder"));
    near_surface.level_set_shape_.writeLevelSet(io_environment);
    InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection, SequencedPolicy> standard_wall_function_correction(water_block_inner, water_block_contact, offset_dist_ref, id_exclude, near_surface);

    //** Build observers in front of the cylinder, 30 cells, 31 bounds *
    for (int i = 0; i < num_observer_points_f + 1; i++)
    {
        monitor_bound_x_f.push_back(x_start_f + i * observe_x_spacing);
    }
    Real monitor_offset = monitor_bound_x_f[num_observer_points_f] - (insert_circle_center[0] - insert_circle_radius);
    for (int i = 0; i < num_observer_points_f + 1; i++)
    {
        monitor_bound_x_f[i] = monitor_bound_x_f[i] - monitor_offset;
    }
    //** Build observers behind the cylinder *
    for (int i = 0; i < num_observer_points_b + 1; i++)
    {
        monitor_bound_x_b.push_back(x_start_b + i * observe_x_spacing);
    }

    SimpleDynamics<fluid_dynamics::GetTimeAverageCenterLineData, SequencedPolicy> get_time_average_center_line_data(water_block_inner, num_observer_points, observe_x_ratio, monitor_bound_y, monitor_bound_x_f, monitor_bound_x_b);

    InteractionDynamics<fluid_dynamics::TurbulentViscousAccelerationWithWall> turbulent_viscous_acceleration(water_block_inner, water_block_contact);
    //InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex_relation);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity> update_eddy_viscosity(water_block);

    /** Turbulent InflowTurbulentCondition.It needs characteristic Length to calculate turbulent length  */
    SimpleDynamics<fluid_dynamics::InflowTurbulentCondition> impose_turbulent_inflow_condition(emitter_buffer, characteristic_length, 0.5);

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> fluid_pressure_force_on_inserted_body(cylinder_contact);
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> fluid_viscous_force_on_inserted_body(cylinder_contact);
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(io_environment, system.real_bodies_);
    ObservedQuantityRecording<Vecd>
        write_fluid_velocity("Velocity", io_environment, fluid_observer_contact);
    RegressionTestTimeAverage<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>>
        write_total_viscous_force_on_inserted_body(io_environment, fluid_viscous_force_on_inserted_body, "TotalViscousForceOnSolid");
    ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>
        write_total_force_on_inserted_body(io_environment, fluid_pressure_force_on_inserted_body, "TotalPressureForceOnSolid");

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    cylinder_normal_direction.exec();
    cylinder.addBodyStateForRecording<Vecd>("NormalDirection");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 200;
    Real output_interval = end_time / 40.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                          /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile();
    get_time_average_center_line_data.output_monitor_x_coordinate();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    int ITER = 0;
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_turbulent_fluid_advection_time_step_size.exec();
            //Real Dt = get_fluid_advection_time_step_size.exec();
            free_stream_surface_indicator.exec();
            update_fluid_density.exec();

            update_eddy_viscosity.exec();
            //viscous_acceleration.exec();
            turbulent_viscous_acceleration.exec();

            transport_velocity_correction.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);

                turbulent_kinetic_energy_acceleration.exec();

                pressure_relaxation.exec(dt);

                emitter_buffer_inflow_condition.exec();

                impose_turbulent_inflow_condition.exec();

                density_relaxation.exec(dt);

                get_velocity_gradient.exec(dt);
                k_equation_relaxation.exec(dt);
                epsilon_equation_relaxation.exec(dt);
                standard_wall_function_correction.exec();

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                //write_body_states.writeToFile();
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
            }
            number_of_iterations++;

            /** inflow injection*/
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
            cylinder_contact.updateConfiguration();

            //if (GlobalStaticVariables::physical_time_ > 0.8)
            //{
            //write_body_states.writeToFile();
            //}

            get_time_average_center_line_data.exec();
            get_time_average_center_line_data.output_time_history_data(end_time * 0.75);
        }

        ITER = ITER + 1;
        //std::cout << "ITER=" << ITER << std::endl;
        //if (GlobalStaticVariables::physical_time_ >=150.)
        //{
        //	output_interval = end_time /  4000.0;
        //}

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_body_states.writeToFile();
        write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
        write_total_force_on_inserted_body.writeToFile(number_of_iterations);
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    get_time_average_center_line_data.get_time_average_data(end_time * 0.75);
    //std::cout << "The time-average data is output " << std::endl;
    if (system.GenerateRegressionData())
    {
        // The lift force at the cylinder is very small and not important in this case.
        write_total_viscous_force_on_inserted_body.generateDataBase({1.0e-2, 1.0e-2}, {1.0e-2, 1.0e-2});
    }
    else
    {
        write_total_viscous_force_on_inserted_body.testResult();
    }
    return 0;
}
