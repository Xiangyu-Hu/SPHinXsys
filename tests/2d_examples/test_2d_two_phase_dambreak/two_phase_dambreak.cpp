/**
 * @file 	two_phase_dambreak.cpp
 * @brief 	2D two-phase dambreak flow.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "two_phase_dambreak.h"
#include "sphinxsys.h"
using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    WaterBlock water_block_shape("WaterBody");
    FluidBody water_block(sph_system, water_block_shape.getName());
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>(water_block_shape);

    AirBlock air_block_shape("AirBody");
    FluidBody air_block(sph_system, air_block_shape.getName());
    air_block.defineMaterial<WeaklyCompressibleFluid>(rho0_a, c_f);
    air_block.generateParticles<BaseParticles, Lattice>(air_block_shape);

    WallBoundary wall_boundary_shape("WallBoundary");
    SolidBody wall_boundary(sph_system, wall_boundary_shape.getName());
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>(wall_boundary_shape);

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_inner(water_block);
    ContactRelation water_air_contact(water_block, {&air_block});
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    InnerRelation air_inner(air_block);
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_wall_contact(air_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, RealBodyVector{&water_block, &air_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_air_complex(water_inner, {&water_air_contact, &water_wall_contact});
    ComplexRelation air_water_complex(air_inner, {&air_water_contact, &air_wall_contact});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SubShapeAndOp *inner_wall_shape_and_op = wall_boundary_shape.getSubShapeAndOpByName("InnerWall");
    SimpleDynamics<NormalDirectionFromSubShapeAndOp> inner_normal_direction(wall_boundary, inner_wall_shape_and_op);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity_to_water(water_block, gravity);
    SimpleDynamics<GravityForce<Gravity>> constant_gravity_to_air(air_block, gravity);
    InteractionDynamics<fluid_dynamics::BoundingFromWall> air_near_wall_bounding(air_wall_contact);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        water_pressure_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        water_density_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        air_pressure_relaxation(air_inner, air_water_contact, air_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        air_density_relaxation(air_inner, air_water_contact, air_wall_contact);

    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface>
        update_water_density_by_summation(water_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_air_density_by_summation(air_inner, air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(air_inner, air_water_contact, air_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_water_advection_time_step_size(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_air_advection_time_step_size(air_block, U_ref);

    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_water_time_step_size(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_air_time_step_size(air_block);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting water_particle_sorting(water_block);
    ParticleSorting air_particle_sorting(air_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection"); // output for debug
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>>
        write_water_mechanical_energy(water_block, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    inner_normal_direction.exec();
    constant_gravity_to_water.exec();
    constant_gravity_to_air.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 10.0;
    Real output_interval = 0.1;
    Real dt = 0.0; /**< Default acoustic time step sizes. */
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
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Force Prior due to viscous force and gravity. */
            time_instance = TickCount::now();

            Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real Dt = SMIN(Dt_f, Dt_a);

            update_water_density_by_summation.exec();
            update_air_density_by_summation.exec();
            air_transport_correction.exec();
            air_near_wall_bounding.exec();

            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
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

                if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
                {
                    write_water_mechanical_energy.writeToFile(number_of_iterations);
                    write_recorded_pressure.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                water_particle_sorting.exec();
                air_particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            air_block.updateCellLinkedList();
            water_air_complex.updateConfiguration();
            air_water_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();

            interval_updating_configuration += TickCount::now() - time_instance;
        }

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
        write_water_mechanical_energy.generateDataBase(1.0e-3);
        write_recorded_pressure.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_water_mechanical_energy.testResult();
        write_recorded_pressure.testResult();
    }

    return 0;
}
