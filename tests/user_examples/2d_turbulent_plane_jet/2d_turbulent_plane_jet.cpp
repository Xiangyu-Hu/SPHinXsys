/**
 * @file 	2d_turbulent_plane_jet.cpp
 * @brief 	2d_turbulent_plane_jet flow with K-Epsilon two equations RANS model.
 * @details This is the one of the basic test cases.
 * @author 	Xiangyu Hu
 */
#include "2d_turbulent_plane_jet.h"

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

    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    //water_block.defineAdaptationRatios(1.4);

    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();

    //wall_boundary.generateParticles<ParticleGeneratorLattice>();
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName())
        : wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(system, "FluidObserver");
    fluid_observer.defineAdaptationRatios(0.0, 1.0);

    for (int j = 0; j < num_observer_points_x; ++j)
    {
        for (int i = 0; i < num_observer_points; ++i)
        {
            observation_locations.push_back(Vecd(x_observe_start + j * observe_spacing_x,
                                                 i * observe_spacing + 0.5 * resolution_ref));
        }
    }
    fluid_observer.generateParticles<ParticleGeneratorObserver>(observation_locations);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ComplexRelation water_block_complex_relation(water_block_inner, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation wall_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_boundary);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(io_environment, {&wall_boundary});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, wall_boundary);
        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(wall_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);

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
        std::cout << "The physics relaxation process of the wall finish !" << std::endl;

        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------

    //Attention! the original one does not use Riemann solver for pressure
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex_relation);
    //Attention! the original one does use Riemann solver for density
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex_relation);

    /** Turbulent.Note: When use wall function, K Epsilon calculation only consider inner */
    InteractionWithUpdate<fluid_dynamics::K_TurbulentModelInner> k_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::GetVelocityGradientInner> get_velocity_gradient(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::E_TurbulentModelInner> epsilon_equation_relaxation(water_block_inner);
    InteractionDynamics<fluid_dynamics::TKEnergyAccComplex> turbulent_kinetic_energy_acceleration(water_block_complex_relation);

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    /** Turbulent standard wall function needs normal vectors of wall. */
    NearShapeSurface near_surface(water_block, makeShared<WallBoundary>("Wall"));
    near_surface.level_set_shape_.writeLevelSet(io_environment);
    InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection, SequencedPolicy> standard_wall_function_correction(water_block_complex_relation, offset_dist_ref, id_exclude, near_surface);

    //SimpleDynamics<fluid_dynamics::GetTimeAverageCrossSectionData,SequencedPolicy> get_time_average_cross_section_data(water_block_inner,num_observer_points);

    InteractionDynamics<fluid_dynamics::TurbulentViscousAccelerationWithWall> turbulent_viscous_acceleration(water_block_complex_relation);
    //InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex_relation);

    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_complex_relation);
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex> inlet_outlet_surface_particle_indicator(water_block_complex_relation);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_complex_relation);
    water_block.addBodyStateForRecording<Real>("Pressure"); // output for debug
    water_block.addBodyStateForRecording<int>("Indicator"); // output for debug

    /** Define the external force for turbulent startup to reduce instability at start-up stage, 1e-4 is from poisulle case */
    SharedPtr<TimeDependentAcceleration> gravity_ptr = makeShared<TimeDependentAcceleration>(Vecd(0.0, 0.0));
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
    //SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);

    /** Turbulent advection time step. */
    ReduceDynamics<fluid_dynamics::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    //ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);

    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity> update_eddy_viscosity(water_block);

    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 50, 0);

    BodyAlignedBoxByCell emitter_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(emitter_buffer);

    /** Turbulent InflowTurbulentCondition.It needs characteristic Length to calculate turbulent length  */
    SimpleDynamics<fluid_dynamics::InflowTurbulentCondition> impose_turbulent_inflow_condition(emitter_buffer, DH, 0.5);

    Vec2d disposer_up_halfsize = Vec2d(0.5 * BW, 0.55 * (DH2 + DE));
    Vec2d disposer_up_translation = Vec2d(DL2 - BW + outlet_length, -0.05 * (DH2 + DE) - DE) + disposer_up_halfsize;
    BodyAlignedBoxByCell disposer_up(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_up_translation)), disposer_up_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_up_outflow_deletion(disposer_up, xAxis);
    ;
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(io_environment, system.real_bodies_);
    //ObservedQuantityRecording<Real> write_fluid_x_velocity("Velocity_X", io_environment, fluid_observer_contact); //For test turbulent model
    //ObservedQuantityRecording<Real> write_fluid_turbu_kinetic_energy("TurbulenceKineticEnergy", io_environment, fluid_observer_contact); //For test turbulent model
    //ObservedQuantityRecording<Real> write_fluid_turbu_dissipation_rate("TurbulentDissipation", io_environment, fluid_observer_contact); //For test turbulent model
    //ObservedQuantityRecording<Real> write_fluid_turbu_viscosity("TurbulentViscosity", io_environment, fluid_observer_contact); //For test turbulent model

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 50.0;
    Real output_interval = end_time / 500.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile();
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
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();

            update_eddy_viscosity.exec();

            //viscous_acceleration.exec();
            turbulent_viscous_acceleration.exec();

            transport_velocity_correction.exec();

            /** Dynamics including pressure relaxation. */
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
                //if(GlobalStaticVariables::physical_time_ >1.83)
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
            disposer_up_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex_relation.updateConfiguration();
            //write_body_states.writeToFile();
            //get_time_average_cross_section_data.exec();
            //get_time_average_cross_section_data.output_cross_section_data();
        }

        ITER = ITER + 1;
        //std::cout << "ITER=" << ITER << std::endl;
        //if (GlobalStaticVariables::physical_time_ >=150.)
        //{
        //	output_interval = end_time /  4000.0;
        //}

        TickCount t2 = TickCount::now();
        write_body_states.writeToFile();

        //write_fluid_x_velocity.writeToFile(); //For test turbulent model
        //write_fluid_turbu_kinetic_energy.writeToFile(); //For test turbulent model
        //write_fluid_turbu_dissipation_rate.writeToFile(); //For test turbulent model
        //write_fluid_turbu_viscosity.writeToFile(); //For test turbulent model

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    write_body_states.writeToFile();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    //get_time_average_cross_section_data.get_time_average_data();
    std::cout << "The time-average data is output " << std::endl;

    return 0;
}
