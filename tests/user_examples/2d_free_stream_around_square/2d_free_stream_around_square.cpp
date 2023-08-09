/**
 * @file 	freestream_flow_around_cylinder.cpp
 * @author 	Yongchuan Yu
 */

#include "sphinxsys.h"
#include "2d_free_stream_around_square.h"
#include "exclusive_shape.h"
#include "fluid_boundary_static_confinement.h"

using namespace SPH;

int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** handle command line arguments. */
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
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
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d::Zero()));
	BodyAlignedBoxByParticle emitter(
		water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
	/** Emitter buffer inflow condition. */
	BodyAlignedBoxByCell emitter_buffer(
		water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
	SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition(emitter_buffer);
	BodyAlignedBoxByCell disposer(
		water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, 0);

	/** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

	/** time-space method to detect surface particles. */
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationInner>
		free_stream_surface_indicator(water_block_inner);

	/** Evaluation of density by freestream approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_fluid_density(water_block_inner);
	/** We can output a method-specific particle data for debug */
	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<int>("SurfaceIndicator");
	
	/** modify the velocity of boundary particles with free-stream velocity. */
	SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint(water_block);
	/** Pressure relaxation. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> pressure_relaxation(water_block_inner);
	/** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
	pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
	/** Density relaxation. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> density_relaxation(water_block_inner);
	/** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationInner> viscous_acceleration(water_block_inner);
	/** Apply transport velocity formulation. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionInner, SequencedPolicy> transport_velocity_correction(water_block_inner);
	/** compute the vorticity. */
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);

    NearShapeSurface near_surface(water_block, makeShared<ExclusiveShape<Cylinder>>("Cylinder"));
    near_surface.level_set_shape_.writeLevelSet(io_environment);
    fluid_dynamics::StaticConfinementGeneral confinement_condition(near_surface);

	free_stream_surface_indicator.post_processes_.push_back(&confinement_condition.free_surface_indication_);
	update_fluid_density.post_processes_.push_back(&confinement_condition.density_summation_);
    pressure_relaxation.post_processes_.push_back(&confinement_condition.pressure_relaxation_);
    density_relaxation.post_processes_.push_back(&confinement_condition.density_relaxation_);
    transport_velocity_correction.post_processes_.push_back(&confinement_condition.transport_velocity_);
    viscous_acceleration.post_processes_.push_back(&confinement_condition.viscous_acceleration_);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	//SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);
	/** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	//InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> fluid_pressure_force_on_inserted_body(cylinder_contact);
	//InteractionDynamics<solid_dynamics::ViscousForceFromFluid> fluid_viscous_force_on_inserted_body(cylinder_contact);
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
	/** initialize configurations for all bodies. */
	sph_system.initializeSystemConfigurations();
	/** computing surface normal direction for the insert body. */
	//cylinder_normal_direction.exec();
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_real_body_states.writeToFile();
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	Real end_time = 50.0;
	Real output_interval = end_time / 400.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TimeInterval interval;
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
			free_stream_surface_indicator.exec();
			update_fluid_density.exec();
			viscous_acceleration.exec();
            transport_velocity_correction.exec(GlobalStaticVariables::physical_time_);

			size_t inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
				/** Fluid pressure relaxation, first half. */
				pressure_relaxation.exec(dt);
				/** Fluid pressure relaxation, second half. */
				density_relaxation.exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				emitter_buffer_inflow_condition.exec();
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
			emitter_inflow_injection.exec();
			disposer_outflow_deletion.exec();

			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_inner.updateConfiguration();
			/** one need update configuration after periodic condition. */
			/** write run-time observation into file */
			//cylinder_contact.updateConfiguration();
		}

		TickCount t2 = TickCount::now();
		/** write run-time observation into file */
		compute_vorticity.exec();
		write_real_body_states.writeToFile();
		//write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
		//write_total_force_on_inserted_body.writeToFile(number_of_iterations);
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
