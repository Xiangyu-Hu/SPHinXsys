/**
 * @file 	2d_turbulent_channel.cpp
 * @brief 	2D_turbulent_channel flow with K-Epsilon two equations RANS model.
 * @details This is the one of the basic test cases.
 * @author 	Xiangyu Hu
 */
//#include "sphinxsys.h"
#include "2d_turbulent_channel_periodicBC.h"

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);
    system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(system, "FluidObserver");
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
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which are only use for now for updating configurations.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------

    //**Attention! the original one does not use Riemann solver for pressure
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    //**Attention! the original one does use Riemann solver for density
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);

    /** Turbulent.Note: When use wall function, K Epsilon and TKE gradient calculation only consider inner */
    //InteractionWithUpdate<fluid_dynamics::K_TurbulentModelInner> k_equation_relaxation(water_block_inner);
    //InteractionWithUpdate<fluid_dynamics::E_TurbulentModelInner> epsilon_equation_relaxation(water_block_inner);
    //InteractionDynamics<fluid_dynamics::TKEnergyAccInner> turbulent_kinetic_energy_acceleration(water_block_inner);

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** Turbulent standard wall function needs normal vectors of wall. */
    //InteractionDynamics<fluid_dynamics::StandardWallFunctionCorrection> standard_wall_function_correction(water_block_complex_relation);

    /** TurbulentViscous cal. includes the molecular one and the eddy viscosity one */
    /** TurbulentViscous cal. uses friction velocity and Y+ that are defined in WallFunction . */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_inner, water_block_contact);
    //InteractionDynamics<fluid_dynamics::TurbulentViscousAccelerationWithWall, SequencedPolicy> turbulent_viscous_acceleration(water_block_complex_relation);

    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);

    //InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
    //	inlet_outlet_surface_particle_indicator(water_block_complex_relation);

    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    water_block.addBodyStateForRecording<Real>("Pressure"); // output for debug
    water_block.addBodyStateForRecording<int>("Indicator"); // output for debug

    /** Define the external force for turbulent startup */
    /**to reduce instability at start-up stage, 1e-4 is from poisulle case */
    //SharedPtr<TimeDependentAcceleration> gravity_ptr = makeShared<TimeDependentAcceleration>(Vecd(1.0e-4, 0.0));
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<Gravity>(Vecd(gravity_g, 0.0)));

    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, water_block.getBodyShapeBounds(), xAxis);

    /** Turbulent advection time step. */
    //ReduceDynamics<fluid_dynamics::TurbulentAdvectionTimeStepSize> get_turbulent_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);

    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /** Turbulent eddy viscosity calculation needs values of Wall Y start. */
    //SimpleDynamics<fluid_dynamics::TurbulentEddyViscosity> update_eddy_viscosity(water_block);

    /**  Try to introduce B correction */
    //InteractionWithUpdate<CorrectedConfigurationInner> correct_configuration(water_block_inner);

    /** Turbulent InflowTurbulentCondition.It needs characteristic Length to calculate turbulent length  */
    //SimpleDynamics<fluid_dynamics::InflowTurbulentCondition> impose_turbulent_inflow_condition(emitter_buffer,DH,0.5);

    //Vec2d disposer_up_halfsize = Vec2d(0.5 * BW, 0.55 * DH);
    //Vec2d disposer_up_translation = Vec2d(DL - BW, -0.05 * DH) + disposer_up_halfsize;
    //BodyAlignedBoxByCell disposer_up(
    //	water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_up_translation)), disposer_up_halfsize));
    //SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_up_outflow_deletion(disposer_up, xAxis);
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

    periodic_condition.update_cell_linked_list_.exec();

    system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 400.0;
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
            //Real Dt = get_turbulent_fluid_advection_time_step_size.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            //inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();

            //update_eddy_viscosity.exec();
            viscous_acceleration.exec();
            //turbulent_viscous_acceleration.exec();

            transport_velocity_correction.exec();
            //** Try to introduce B correction *
            //correct_configuration.exec();
            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);

                //turbulent_kinetic_energy_acceleration.exec();

                pressure_relaxation.exec(dt);

                //impose_turbulent_inflow_condition.exec();

                density_relaxation.exec(dt);

                //k_equation_relaxation.exec(dt);
                //epsilon_equation_relaxation.exec(dt);
                //standard_wall_function_correction.exec();

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

            /** Periodic condition. */
            periodic_condition.bounding_.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);

            /** Periodic condition. */
            periodic_condition.update_cell_linked_list_.exec();

            water_block_complex.updateConfiguration();
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

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    return 0;
}
