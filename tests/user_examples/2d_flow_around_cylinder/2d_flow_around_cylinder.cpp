/**
 * @file 	2d_flow_around_cylinder.cpp
 * @brief 	This is the benchmark test for the wall modeling of viscous flow.
 * @details We consider a flow passing by a cylinder in 2D.
 * @author 	Yongchuan Yu
 */
#include "2d_flow_around_cylinder.h"
#include "sphinxsys.h"
#include "level_set_confinement.h"

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
// handle command line arguments
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);
    ParameterizationIO &parameterization_io = io_environment.defineParameterizationIO();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineParticlesAndMaterial<BaseParticles, ParameterizedWaterMaterial>(parameterization_io, rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_locations);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Initialize particle acceleration. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, water_block.getBodyShapeBounds(), xAxis);
    /** Periodic BCs in y direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_block, water_block.getBodyShapeBounds(), yAxis);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationInner> update_density_by_summation(water_block_inner);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation using Verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> pressure_relaxation(water_block_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalf> density_relaxation(water_block_inner);
    /** Computing viscous acceleration with wall. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationInner> viscous_acceleration(water_block_inner);
    /** Impose transport velocity. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionInner> transport_velocity_correction(water_block_inner);
    /** Computing vorticity in the flow. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    /** free stream boundary condition. */
    BodyRegionByCell free_stream_buffer(water_block, makeShared<MultiPolygonShape>(createBufferShape()));
    SimpleDynamics<FreeStreamCondition> freestream_condition(free_stream_buffer);
    /** time-space method to detect surface particles. */
  /*  InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationInner>
		free_stream_surface_indicator(water_block_inner);*/

    NearShapeSurface near_surface(water_block, makeShared<InverseShape<Cylinder>>("Cylinder"));
    near_surface.level_set_shape_.writeLevelSet(io_environment);
    fluid_dynamics::StaticConfinementGeneral confinement_condition(near_surface);

    //free_stream_surface_indicator.post_processes_.push_back(&confinement_condition.free_surface_indication_);
	update_density_by_summation.post_processes_.push_back(&confinement_condition.density_summation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition.density_relaxation_);
    transport_velocity_correction.post_processes_.push_back(&confinement_condition.transport_velocity_);
    viscous_acceleration.post_processes_.push_back(&confinement_condition.viscous_acceleration_);

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
   /* InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_cylinder(cylinder_contact);
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> pressure_force_on_cylinder(cylinder_contact);*/
    /** Computing viscous force acting on wall with wall model. */
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    
    ObservedQuantityRecording<Vecd>
        write_fluid_velocity("Velocity", io_environment, fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition_x.update_cell_linked_list_.exec();
    periodic_condition_y.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
   
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 200.0;
    Real output_interval = end_time / 200.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;

        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            //free_stream_surface_indicator.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec(GlobalStaticVariables::physical_time_);

            /** FSI for viscous force. */
            
            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation, first half. */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
       
                /** Fluid pressure relaxation, second half. */
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                freestream_condition.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition_x.bounding_.exec();
            periodic_condition_y.bounding_.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            periodic_condition_x.update_cell_linked_list_.exec();
            periodic_condition_y.update_cell_linked_list_.exec();
            /** one need update configuration after periodic condition. */
            water_block_inner.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;


    return 0;
}