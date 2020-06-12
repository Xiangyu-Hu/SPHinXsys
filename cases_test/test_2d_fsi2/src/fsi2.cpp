/**
 * @file 	fsi2.cpp
 * @brief 	This is the benchmark test of fliud-structure interaction.
 * @details We consider a flow-induced vibration of an elastic beam behind a cylinder in 2D.
 *			The case can be found in Chi Zhang, Massoud Rezavand, Xiangyu Hu,
 *			Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 *			Journal of Computation Physics 404 (2020) 109135.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 * @version 0.1
 */
#include "sphinxsys.h"

/** case file to setup the test case */
#include "fsi2_case.h"

using namespace SPH;

int main(int ac, char* av[])
{
	/** Build up -- a SPHSystem -- */
	SPHSystem system(Vec2d(-DLsponge - BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;

	//handle command line arguments
	system.handleCommandlineOptions(ac, av);

	/**
	 * @brief Creating body, materials and particles for a water block.
	 */
	WaterBlock* water_block
		= new WaterBlock(system, "WaterBody", 0, ParticlesGeneratorOps::lattice);
	WaterMaterial* water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Creating body and particles for the wall boundary.
	 */
	WallBoundary* wall_boundary
		= new WallBoundary(system, "Wall", 0, ParticlesGeneratorOps::lattice);
	SolidParticles 	solid_particles(wall_boundary);
	/**
	 * @brief 	Creating body, materials and particles for the elastic beam (inserted body).
	 */
	InsertedBody* inserted_body = new InsertedBody(system, "InsertedBody", 1, ParticlesGeneratorOps::lattice);
	InsertBodyMaterial* insert_body_material = new InsertBodyMaterial();
	ElasticSolidParticles 	inserted_body_particles(inserted_body, insert_body_material);
	/**
	 * @brief 	Particle and body creation of beam and fluid observers.
	 */
	BeamObserver* beam_observer = new BeamObserver(system, "BeamObserver", 0, ParticlesGeneratorOps::direct);
	BaseParticles 				beam_observer_particles(beam_observer);
	FluidObserver* fluid_observer = new FluidObserver(system, "FluidObserver", 0, ParticlesGeneratorOps::direct);
	BaseParticles				flow_observer_particles(fluid_observer);
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	In_Output 							in_output(system);
	//WriteBodyStatesToVtu 				write_real_body_states(in_output, system.real_bodies_);
	WriteBodyStatesToPlt 				write_real_body_states(in_output, system.real_bodies_);
	WriteRestart						write_restart_files(in_output, system.real_bodies_);
	ReadRestart							read_restart_files(in_output, system.real_bodies_);
	/**
	 * @brief 	Body contact map.
	 * @details The contact map gives the topological conntections between the bodies.
	 * 			Basically the the range of bodies to build neighbor particle lists.
	 */
	SPHBodyTopology body_topology = { { water_block, { wall_boundary, inserted_body } },
									  { wall_boundary, { } }, { inserted_body, { water_block } },
									  { beam_observer, {inserted_body} }, { fluid_observer, { water_block } } };
	system.SetBodyTopology(&body_topology);
	/**
	 * @brief 	Methods used for general methods.
	 */
	 /** Update the cell linked list of a body. */
	ParticleDynamicsCellLinkedList 			update_water_block_cell_linked_list(water_block);
	/** Update the cell linked list of a body. */
	ParticleDynamicsCellLinkedList 			update_inserted_body_cell_linked_list(inserted_body);
	/** Update all configurations (including inner and contact) of a body. */
	ParticleDynamicsConfiguration 			update_water_block_configuration(water_block);
	/** Update the contact configuration for a body. */
	ParticleDynamicsContactConfiguration 	update_inserted_body_contact_configuration(inserted_body);
	/** Update the contact configuration for the flow observer. */
	ParticleDynamicsContactConfiguration 	update_fluid_observer_body_contact_configuration(fluid_observer);
	/** Update inner configuration of a body. */
	ParticleDynamicsInnerConfiguration update_inserted_body_inner_configuration(inserted_body);
	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(water_block, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(water_block, 0);

	/** check whether run particle relaxation for body fiited particle distribution. */
	if (system.run_particle_relaxation_) 
	{
		/** add background level set for particle realxation. */
		inserted_body->addBackgroundMesh();
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		 /** Random reset the insert body particle position. */
		RandomizePartilePosition  random_inserted_body_particles(inserted_body);
		/** Write backgroung level set. */
		WriteBodyMeshToPlt write_inserted_body_background_mesh(in_output, inserted_body);
		/** Write the body state to Vtu file. */
		WriteBodyStatesToVtu 		write_inserted_body_to_vtu(in_output, { inserted_body });
		/** Write the particle reload files. */
		WriteReloadParticle 		write_particle_reload_files(in_output, { inserted_body });

		/** bounding particles to insert body surface. */
		relax_dynamics::BodySurfaceBounding
			body_surface_bounding(inserted_body, new NearBodySurface(inserted_body));
		/** Compute the time step for particle relaxation. */
		relax_dynamics::GetTimeStepSize get_solid_relax_timestep(inserted_body);
		/** Physics relaxation algorithm. */
		relax_dynamics::PhysicsRelaxationInner 	relax_process_for_solid(inserted_body);
		/** finilaizing  particle number desnity and inital position after relaxatoin. */
		relax_dynamics::FinalizingParticleRelaxation finalizing_inserted_body_particles(inserted_body);
		/**
		  * @brief 	Particle relaxation starts here.
		  */
		write_inserted_body_background_mesh.WriteToFile(0.0);
		update_inserted_body_cell_linked_list.parallel_exec();
		update_inserted_body_inner_configuration.parallel_exec();
		body_surface_bounding.parallel_exec();
		random_inserted_body_particles.parallel_exec(0.25);
		write_real_body_states.WriteToFile(0.0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		Real dt_p = 0.0;
		while (ite_p < 1000)
		{
			dt_p = get_solid_relax_timestep.parallel_exec();

			relax_process_for_solid.parallel_exec(dt_p);
			body_surface_bounding.parallel_exec();
			ite_p += 1;

			update_inserted_body_cell_linked_list.parallel_exec();
			update_inserted_body_inner_configuration.parallel_exec();
			if (ite_p % 200 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_inserted_body_to_vtu.WriteToFile(Real(ite_p) * 1.0e-4);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
		finalizing_inserted_body_particles.parallel_exec();

		/** Output results. */
		write_particle_reload_files.WriteToFile(0.0);
		return 0;
	}
	/**
	 * @brief 	Define all numerical methods which are used in FSI.
	 */
	 /** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_boundary, {});
	/** Initialize normal direction of the inserted body. */
	solid_dynamics::NormalDirectionSummation 	get_inserted_body_normal(inserted_body, {});
	/** Corrected strong configuration for the elastic  insertbody. */
	solid_dynamics::CorrectConfiguration 		inserted_body_corrected_configuration_in_strong_form(inserted_body);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_desnity(water_block, { wall_boundary, inserted_body });
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize 	get_fluid_adevction_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationFirstHalf
		pressure_relaxation_first_half(water_block, { wall_boundary, inserted_body });
	fluid_dynamics::PressureRelaxationSecondHalf
		pressure_relaxation_second_half(water_block, { wall_boundary, inserted_body });
	/** Computing viscous acceleration. */
	fluid_dynamics::ComputingViscousAcceleration 	viscous_acceleration(water_block, { wall_boundary, inserted_body });
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrection 	transport_velocity_correction(water_block, { wall_boundary, inserted_body });
	/** Computing vorticity in the flow. */
	fluid_dynamics::ComputingVorticityInFluidField 	compute_vorticity(water_block);
	/** Inflow boundary condition. */
	ParabolicInflow		parabolic_inflow(water_block, new InflowBuffer(water_block, "Buffer"));
	/**
	 * @brief Algorithms of FSI.
	 */
	 /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidPressureForceOnSolid 	fluid_pressure_force_on_insrted_body(inserted_body, { water_block });
	solid_dynamics::FluidViscousForceOnSolid 	fluid_viscous_force_on_insrted_body(inserted_body, { water_block });
	/** Computing the average velocity. */
	solid_dynamics::InitializeDisplacement 			inserted_body_initialize_displacement(inserted_body);
	solid_dynamics::UpdateAverageVelocity 			inserted_body_average_velocity(inserted_body);
	/**
	 * @brief Algorithms of solid dynamics.
	 */
	 /** Compute time step size of elastic solid. */
	solid_dynamics::GetAcousticTimeStepSize 	inserted_body_computing_time_step_size(inserted_body);
	/** Stress relaxation for the inserted body. */
	solid_dynamics::StressRelaxationFirstHalf 	inserted_body_stress_relaxation_first_half(inserted_body);
	solid_dynamics::StressRelaxationSecondHalf 	inserted_body_stress_relaxation_second_half(inserted_body);
	/** Constrain region of the inserted body. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_beam_base(inserted_body, new BeamBase(inserted_body, "BeamBase"));
	/** Update norm .*/
	solid_dynamics::UpdateElasticNormalDirection 	inserted_body_update_normal(inserted_body);
	/**
	 * @brief Write observation data into files.
	 */
	WriteTotalViscousForceOnSolid 		write_total_viscous_force_on_inserted_body(in_output, inserted_body);

	WriteAnObservedQuantity<Vecd, BaseParticles,
		BaseParticleData, &BaseParticles::base_particle_data_, &BaseParticleData::pos_n_>
		write_beam_tip_displacement("Displacement", in_output, beam_observer, inserted_body);
	WriteAnObservedQuantity<Vecd, BaseParticles,
		BaseParticleData, &BaseParticles::base_particle_data_, &BaseParticleData::vel_n_>
		write_fluid_velocity("Velocity", in_output, fluid_observer, water_block);

	/**
	 * @brief Pre-simulation.
	 */
	 /** Using relaxed particle distribution if needed. */
	if (system.reload_particles_) {
		unique_ptr<ReadReloadParticle>	
			reload_insert_body_particles(new ReadReloadParticle(in_output, { inserted_body }, { "InsertedBody" }));
		reload_insert_body_particles->ReadFromFile();
	}
	/** intialize cell linked lists for all bodies. */
	system.InitializeSystemCellLinkedLists();
	/** peroidic condition applied after the mesh cell linked list build up
	  * but before the configuration bulid up. */
	periodic_condition.parallel_exec();
	/** intialize configurations for all bodies. */
	system.InitializeSystemConfigurations();
	/** computing surface normal direction for the wall. */
	get_wall_normal.parallel_exec();
	/** computing surface normal direction for the insert body. */
	get_inserted_body_normal.parallel_exec();
	/** computing linear reproducing configuration for the insert body. */
	inserted_body_corrected_configuration_in_strong_form.parallel_exec();

	/**
	 * @brief The time stepping starts here.
	 */
	if (system.restart_step_ != 0) {
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(system.restart_step_);
		update_inserted_body_cell_linked_list.parallel_exec();
		update_water_block_cell_linked_list.parallel_exec();
		periodic_condition.parallel_exec();
		/** one need update configuration after peroidic condition. */
		update_water_block_configuration.parallel_exec();
		update_inserted_body_contact_configuration.parallel_exec();
		inserted_body_update_normal.parallel_exec();
	}
	/** first output*/
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

	int number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 200.0;			/**< End time. */
	Real D_Time = End_Time / 200.0;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
	size_t inner_ite_dt = 0;
	size_t inner_ite_dt_s = 0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * @brief Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integeral_time < D_Time) {
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_viscous_force_on_insrted_body.parallel_exec();
			/** Update normal direction on elastic body.*/
			inserted_body_update_normal.parallel_exec();
			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				/** Fluid pressure relaxation, first half. */
				pressure_relaxation_first_half.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_insrted_body.parallel_exec();
				/** Fluid pressure relaxation, second half. */
				pressure_relaxation_second_half.parallel_exec(dt);

				/** Solid dynamics. */
				inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				inserted_body_initialize_displacement.parallel_exec();
				while (dt_s_sum < dt) {

					dt_s = inserted_body_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;
					inserted_body_stress_relaxation_first_half.parallel_exec(dt_s);
					constrain_beam_base.parallel_exec();
					inserted_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
					inner_ite_dt_s++;
				}
				inserted_body_average_velocity.parallel_exec(dt);

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					write_restart_files.WriteToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_bounding.parallel_exec();
			update_water_block_cell_linked_list.parallel_exec();
			periodic_condition.parallel_exec();
			update_water_block_configuration.parallel_exec();
			/** Inserted body contact configuration. */
			update_inserted_body_cell_linked_list.parallel_exec();
			update_inserted_body_contact_configuration.parallel_exec();
			/** write run-time observation into file */
			write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		write_total_viscous_force_on_inserted_body.WriteToFile(GlobalStaticVariables::physical_time_);
		update_fluid_observer_body_contact_configuration.parallel_exec();
		write_fluid_velocity.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}

