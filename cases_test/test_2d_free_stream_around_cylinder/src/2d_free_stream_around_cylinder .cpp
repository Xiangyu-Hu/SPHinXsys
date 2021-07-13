/**
* @file 	freestream_flow_around_cylinder_case.h
* @brief 	This is the case file for the test of free-stream flow.
* @details  We consider a flow pass the cylinder with freestream boundary condition in 2D.
* @author 	Xiangyu Hu, Shuoguo Zhang
*/

#include "sphinxsys.h"

/** case file to setup the test case */
#include "2d_free_stream_around_cylinder.h"

using namespace SPH;

int main(int ac, char* av[])
{
	/** Build up -- a SPHSystem -- */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);
	/** Define the external force. */
	TimeDependentAcceleration gravity(Vec2d(0.0, 0.0));
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;

	/** handle command line arguments. */
	system.handleCommandlineOptions(ac, av);
	In_Output in_output(system);

	/**
	 * @brief Creating body, materials and particles for a water block.
	 */
	WaterBlock* water_block = new WaterBlock(system, "WaterBody");
	water_block->setBodyDomainBounds(fluid_body_domain_bounds);
	WaterMaterial* water_material = new WaterMaterial();
	FluidParticles fluid_particles(water_block, water_material);
	/**
	 * @brief 	Creating body, materials and particles for solid cylinder.
	 */
	Cylinder* cylinder = new Cylinder(system, "Cylinder");
	if (!system.run_particle_relaxation_ && system.reload_particles_) cylinder->useParticleGeneratorReload();
	SolidParticles cylinder_particles(cylinder);
	/**
	 * @brief 	Particle and body creation of beam and fluid observers.
	 */
	FluidObserver* fluid_observer = new FluidObserver(system, "FluidObserver");
	BaseParticles flow_observer_particles(fluid_observer);
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	BodyStatesRecordingToVtu write_real_body_states(in_output, system.real_bodies_);
	RestartIO restart_io(in_output, system.real_bodies_);
	/** topology */
	BodyRelationInner* water_block_inner = new BodyRelationInner(water_block);
	ComplexBodyRelation* water_block_complex = new ComplexBodyRelation(water_block_inner, { cylinder });
	BodyRelationContact* cylinder_contact = new BodyRelationContact(cylinder, { water_block });
	BodyRelationContact* fluid_observer_contact = new BodyRelationContact(fluid_observer, { water_block });

	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_)
	{
		/** body topology only for particle relaxation */
		BodyRelationInner* cylinder_inner = new BodyRelationInner(cylinder);
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		 /** Random reset the insert body particle position. */
		RandomizePartilePosition random_inserted_body_particles(cylinder);
		/** Write the body state to Vtu file. */
		BodyStatesRecordingToVtu write_inserted_body_to_vtu(in_output, { cylinder });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, { cylinder });

		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner);
		/**
		  * @brief 	Particle relaxation starts here.
		  */
		random_inserted_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_real_body_states.writeToFile(0.0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_inserted_body_to_vtu.writeToFile(Real(ite_p) * 1.0e-4);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.writeToFile(0.0);
		return 0;
	}

	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block, &gravity);
	/** Emitter condition. */
	EmitterBuffer* emitter = new EmitterBuffer(water_block, "Emitter");
	EmitterInflowCondition inlet_outlet_inflow(water_block, emitter);
	/** Injection condition. */
	InjectionBuffer* injection = new InjectionBuffer(water_block, "Injection");
	fluid_dynamics::EmitterInflowInjecting inflow_emitter(water_block, injection, 50, 0, true);
	/** time-space method to detect surface particles. */
	fluid_dynamics::SurfaceParticlesIndicator free_stream_surface_indicator(water_block_complex,water_block_inner,0);
	/** Evaluation of density by freestream approach. */
	fluid_dynamics::DensitySummationFreeStreamComplex update_fluid_density(water_block_complex);
	/** We can output a method-specific particle data for debug reason */
	fluid_particles.addAVariableToWrite<indexScalar, Real>("Pressure");
	fluid_particles.addAVariableToWrite<indexInteger, Real>("SurfaceIndicator");
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** modify the velocity of boundary particles with far-field velocity. */
	fluid_dynamics::FreeStreamBoundaryVelocityCorrection velocity_boundary_condition_constraint(water_block_inner);
	/** Pressure relaxation. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex);
	/** correct the velotity of boundary particles with far-field velocity through the post process of pressure relaxation. */
	pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
	/** Density relaxation. */
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** recycle real fluid particle to buffer particles at outlet. */
	OpenBoundaryConditionInAxisDirection tansfer_to_buffer_particles_upper_bound(water_block, xAxis, positiveDirection);
	/** compute the vorticity. */
	fluid_dynamics::VorticityInner compute_vorticity(water_block_inner);
	/**
	 * @brief Algorithms of FSI.
	 */
	 /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_inserted_body(cylinder_contact);
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_inserted_body(cylinder_contact);
	/**
	 * @brief Write observation data into files.
	 */
	ObservedQuantityRecording<indexVector, Vecd>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);
	BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid> write_total_viscous_force_on_inserted_body(in_output, cylinder);
	BodyReducedQuantityRecording<solid_dynamics::TotalForceOnSolid> write_total_force_on_inserted_body(in_output, cylinder);
	/**
	 * @brief Pre-simulation.
	 */
	 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	  * but before the configuration build up. */
	  /** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the insert body. */
	cylinder_particles.initializeNormalDirectionFromGeometry();
	free_stream_surface_indicator.first_layer.parallel_exec();
	free_stream_surface_indicator.second_and_third_layers.parallel_exec();
	free_stream_surface_indicator.inlet_outlet_layers.parallel_exec();

	/**
	 * @brief The time stepping starts here.
	 */
	if (system.restart_step_ != 0) {
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		cylinder->updateCellLinkedList();
		water_block->updateCellLinkedList();
		/** one need update configuration after periodic condition. */
		water_block_complex->updateConfiguration();
		cylinder_contact->updateConfiguration();
	}
	/** first output*/
	write_real_body_states.writeToFile(GlobalStaticVariables::physical_time_);

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 200.0;			/**< End time. */
	Real D_Time = End_Time / 400.0;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */
	size_t inner_ite_dt = 0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * @brief Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_viscous_force_on_inserted_body.parallel_exec();
			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				/** Fluid pressure relaxation, first half. */
				pressure_relaxation.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_inserted_body.parallel_exec();
				/** Fluid pressure relaxation, second half. */
				density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				inlet_outlet_inflow.parallel_exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			inflow_emitter.exec();
			tansfer_to_buffer_particles_upper_bound.particle_type_transfer.parallel_exec();

			water_block->updateCellLinkedList();
			water_block_complex->updateConfiguration();
			/** one need update configuration after periodic condition. */
			/** write run-time observation into file */
			cylinder_contact->updateConfiguration();
			cylinder_particles.initializeNormalDirectionFromGeometry();
			free_stream_surface_indicator.first_layer.parallel_exec();
			free_stream_surface_indicator.second_and_third_layers.parallel_exec();
			free_stream_surface_indicator.inlet_outlet_layers.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.writeToFile();
		write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
		write_total_force_on_inserted_body.writeToFile(number_of_iterations);
		fluid_observer_contact->updateConfiguration();
		write_fluid_velocity.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}

