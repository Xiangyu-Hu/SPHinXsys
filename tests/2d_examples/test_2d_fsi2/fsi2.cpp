/**
 * @file 	fsi2.cpp
 * @brief 	This is the benchmark test of fluid-structure interaction.
 * @details We consider a flow-induced vibration of an elastic beam behind a cylinder in 2D.
 *			The case can be found in Chi Zhang, Massoud Rezavand, Xiangyu Hu,
 *			Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 *			Journal of Computation Physics 404 (2020) 109135.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */
#include "sphinxsys.h"

#include "fsi2.h" //	case file to setup the test case
using namespace SPH;

int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;
//handle command line arguments
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	/** output environment. */
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	WaterBlock water_block(system, "WaterBody");
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));

	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);

	InsertedBody inserted_body(system, "InsertedBody");
	SharedPtr<ParticleGenerator> inserted_body_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		inserted_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, inserted_body.getBodyName());
	ElasticSolidParticles inserted_body_particles(inserted_body, 
		makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson), inserted_body_particle_generator);

	ObserverBody beam_observer(system, "BeamObserver");
	ObserverParticles beam_observer_particles(beam_observer, makeShared<BeamObserverParticleGenerator>());

	ObserverBody fluid_observer(system, "FluidObserver");
	ObserverParticles flow_observer_particles(fluid_observer, makeShared<FluidObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner inserted_body_inner(inserted_body);
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary, &inserted_body});
	BodyRelationContact inserted_body_contact(inserted_body, {&water_block});
	BodyRelationContact beam_observer_contact(beam_observer, {&inserted_body});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		RandomizePartilePosition random_inserted_body_particles(inserted_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_inserted_body_to_vtp(in_output, {&inserted_body});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, {&inserted_body});
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(inserted_body_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_inserted_body_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
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
	/** Corrected configuration for the elastic  insertbody. */
	solid_dynamics::CorrectConfiguration inserted_body_corrected_configuration(inserted_body_inner);
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationComplex update_density_by_summation(water_block_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInner compute_vorticity(water_block_complex.inner_relation_);
	/** Inflow boundary condition. */
	MultiPolygonShape inflow_buffer_shape(createInflowBufferShape());
	BodyRegionByCell inflow_buffer(water_block, "Buffer", inflow_buffer_shape);
	ParabolicInflow parabolic_inflow(water_block, inflow_buffer);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition(water_block, xAxis);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	/** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidForceOnSolidUpdate fluid_force_on_solid_update(inserted_body_contact);
	/** Compute the average velocity of the insert body. */
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(inserted_body);
	//----------------------------------------------------------------------
	//	Algorithms of solid dynamics.
	//----------------------------------------------------------------------
	/** Compute time step size of elastic solid. */
	solid_dynamics::AcousticTimeStepSize inserted_body_computing_time_step_size(inserted_body);
	/** Stress relaxation for the inserted body. */
	solid_dynamics::StressRelaxationFirstHalf inserted_body_stress_relaxation_first_half(inserted_body_inner);
	solid_dynamics::StressRelaxationSecondHalf inserted_body_stress_relaxation_second_half(inserted_body_inner);
	/** Constrain region of the inserted body. */
	MultiPolygonShape beam_base_shape(createBeamBaseShape());
	BodyRegionByParticle beam_base(inserted_body, "BeamBase", beam_base_shape);
	solid_dynamics::ConstrainSolidBodyRegion constrain_beam_base(inserted_body, beam_base);
	/** Update norm .*/
	solid_dynamics::UpdateElasticNormalDirection inserted_body_update_normal(inserted_body);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states(in_output, system.real_bodies_);
	RestartIO restart_io(in_output, system.real_bodies_);
	RegressionTestTimeAveraged<BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid>>
		write_total_viscous_force_on_inserted_body(in_output, inserted_body);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_beam_tip_displacement("Position", in_output, beam_observer_contact);
	ObservedQuantityRecording<Vecd>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	/** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	  * but before the configuration build up. */
	periodic_condition.update_cell_linked_list_.parallel_exec();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the wall. */
	wall_particles.initializeNormalDirectionFromBodyShape();
	/** computing surface normal direction for the insert body. */
	inserted_body_particles.initializeNormalDirectionFromBodyShape();
	/** computing linear reproducing configuration for the insert body. */
	inserted_body_corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		inserted_body.updateCellLinkedList();
		water_block.updateCellLinkedList();
		periodic_condition.update_cell_linked_list_.parallel_exec();
		/** one need update configuration after periodic condition. */
		water_block_complex.updateConfiguration();
		inserted_body_contact.updateConfiguration();
		inserted_body_update_normal.parallel_exec();
	}
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 200.0;			/**< End time. */
	Real D_Time = End_Time / 200.0; /**< time stamps for output. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_real_body_states.writeToFile();
	write_beam_tip_displacement.writeToFile(number_of_iterations);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_force_on_solid_update.viscous_force_.parallel_exec();
			/** Update normal direction on elastic body.*/
			inserted_body_update_normal.parallel_exec();
			size_t inner_ite_dt = 0;
			size_t inner_ite_dt_s = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt);
				/** Fluid pressure relaxation */
				pressure_relaxation.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_force_on_solid_update.parallel_exec();
				/** Fluid density relaxation */
				density_relaxation.parallel_exec(dt);

				/** Solid dynamics. */
				inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					Real dt_s = SMIN(inserted_body_computing_time_step_size.parallel_exec(), dt - dt_s_sum);
					inserted_body_stress_relaxation_first_half.parallel_exec(dt_s);
					constrain_beam_base.parallel_exec();
					inserted_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
					inner_ite_dt_s++;
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_condition.bounding_.parallel_exec();

			water_block.updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
			water_block_complex.updateConfiguration();
			/** one need update configuration after periodic condition. */
			inserted_body.updateCellLinkedList();
			inserted_body_contact.updateConfiguration();
			/** write run-time observation into file */
			write_beam_tip_displacement.writeToFile(number_of_iterations);
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.writeToFile();
		write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
		fluid_observer_contact.updateConfiguration();
		write_fluid_velocity.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_total_viscous_force_on_inserted_body.newResultTest();
	write_beam_tip_displacement.newResultTest();

	return 0;
}
