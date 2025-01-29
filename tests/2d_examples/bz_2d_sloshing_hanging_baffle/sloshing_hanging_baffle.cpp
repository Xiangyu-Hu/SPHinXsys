/**
 * @file sloshing_hanging_baffle.cpp
 * @author Bo Zhang and Xiangyu Hu
 */
#include "sloshing_hanging_baffle.h" // case file to setup the test case
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody baffle(sph_system, makeShared<Baffle>("Baffle"));
    baffle.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    baffle.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of gate observer.
    //----------------------------------------------------------------------
    ObserverBody baffle_observer(sph_system, "Observer");
    baffle_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation baffle_inner(baffle);
    ContactRelation water_block_contact(water_block, { &wall_boundary, &baffle });
    ContactRelation baffle_contact(baffle, { &water_block });
    ContactRelation baffle_observer_contact(baffle_observer, { &baffle });
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    // For typical fluid-structure interaction, we first define structure dynamics,
    // Then fluid dynamics and the corresponding coupling dynamics.
    // The coupling with multi-body dynamics will be introduced at last.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> baffle_normal_direction(baffle);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> baffle_corrected_configuration(baffle_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> baffle_stress_relaxation_first_half(baffle_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> baffle_stress_relaxation_second_half(baffle_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> baffle_computing_time_step_size(baffle);
    BodyRegionByParticle baffle_constraint_part(baffle, makeShared<MultiPolygonShape>(BaffleConstraint()));
    SimpleDynamics<FixBodyPartConstraint> baffle_constraint(baffle_constraint_part);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    VariableGravity variablegravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<VariableGravity>> variable_gravity(water_block, variablegravity);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(baffle);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(baffle);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(baffle_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_baffle_tip_displacement("Position", baffle_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    wall_boundary_normal_direction.exec();
    /** computing surface normal direction for the insert body. */
    baffle_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    baffle_corrected_configuration.exec();
    variable_gravity.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    write_baffle_tip_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 11.0;
    Real output_interval = end_time / 50.0;
    Real dt = 0.0;   /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    Real total_time = 0.0;
    Real relax_time = 1.0;
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            variable_gravity.exec();
            update_density_by_summation.exec();
            /** Update correction matrix for fluid */
            //corrected_configuration_fluid.exec();

            /** Update normal direction on elastic body.*/
            insert_body_update_normal.exec();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    baffle_stress_relaxation_first_half.exec(dt_s);
                    baffle_constraint.exec();
                    baffle_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    dt_s = baffle_computing_time_step_size.exec();
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                total_time += dt;
                if(total_time >= relax_time)
                    physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << " Total Time = " << total_time
                          << " Physical Time = " << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** Water block configuration. */
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            /** one need update configuration after periodic condition. */
            baffle.updateCellLinkedList();
            baffle_contact.updateConfiguration();
            if (total_time >= relax_time)
            {
                /** write run-time observation into file */
                write_baffle_tip_displacement.writeToFile(number_of_iterations);
            }
            
        }
        TickCount t2 = TickCount::now();
        if (total_time >= relax_time)
        {
            /** write run-time observation into file */
            //write_baffle_tip_displacement.writeToFile(number_of_iterations);
            write_real_body_states.writeToFile();
        }
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
