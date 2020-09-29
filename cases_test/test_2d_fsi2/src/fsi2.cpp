/**
 * @file 	fsi2.cpp
 * @brief 	This is the benchmark test of fluid-structure interaction.
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
	SPHSystem system(Vec2d(- DL_sponge - BW, - BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;

	//handle command line arguments
	system.handleCommandlineOptions(ac, av);

	/**
	 * @brief Creating body, materials and particles for a water block.
	 */
	WaterBlock* water_block	= new WaterBlock(system, "WaterBody", 0);
	WaterMaterial* water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Creating body and particles for the wall boundary.
	 */
	WallBoundary* wall_boundary	= new WallBoundary(system, "Wall", 0);
	SolidParticles 	solid_particles(wall_boundary);
	/**
	 * @brief 	Creating body, materials and particles for the elastic beam (inserted body).
	 */
	InsertedBody* inserted_body = new InsertedBody(system, "InsertedBody", 1);
	InsertBodyMaterial* insert_body_material = new InsertBodyMaterial();
	ElasticSolidParticles 	inserted_body_particles(inserted_body, insert_body_material);
	/**
	 * @brief 	Particle and body creation of beam and fluid observers.
	 */
	BeamObserver* beam_observer = new BeamObserver(system, "BeamObserver", 0);
	BaseParticles 				beam_observer_particles(beam_observer);
	FluidObserver* fluid_observer = new FluidObserver(system, "FluidObserver", 0);
	BaseParticles				flow_observer_particles(fluid_observer);
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToVtu 				write_real_body_states(in_output, system.real_bodies_);
	WriteRestart						write_restart_files(in_output, system.real_bodies_);
	ReadRestart							read_restart_files(in_output, system.real_bodies_);
	/** topology */
	SPHBodyInnerRelation* water_block_inner = new SPHBodyInnerRelation(water_block);
	SPHBodyInnerRelation* inserted_body_inner = new SPHBodyInnerRelation(inserted_body);
	SPHBodyComplexRelation* water_block_complex = new SPHBodyComplexRelation(water_block_inner, { wall_boundary, inserted_body });
	SPHBodyComplexRelation* wall_complex = new SPHBodyComplexRelation(wall_boundary, {});
	SPHBodyComplexRelation* inserted_body_complex = new SPHBodyComplexRelation(inserted_body_inner, {});
	SPHBodyContactRelation* inserted_body_contact = new SPHBodyContactRelation(inserted_body, { water_block });
	SPHBodyContactRelation* beam_observer_contact = new SPHBodyContactRelation(beam_observer, { inserted_body });
	SPHBodyContactRelation* fluid_observer_contact = new SPHBodyContactRelation(fluid_observer, { water_block });
	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(water_block, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(water_block, 0);

	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_) 
	{
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		/** Random reset the insert body particle position. */
		RandomizePartilePosition  random_inserted_body_particles(inserted_body);
		/** Write the body state to Vtu file. */
		WriteBodyStatesToVtu 		write_inserted_body_to_vtu(in_output, { inserted_body });
		/** Write the particle reload files. */
		WriteReloadParticle 		write_particle_reload_files(in_output, { inserted_body });

		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(inserted_body_inner);
		/** finalizing  particle number density and inital position after relaxatoin. */
		relax_dynamics::FinalizingParticleRelaxation finalizing_inserted_body_particles(inserted_body);
		/**
		  * @brief 	Particle relaxation starts here.
		  */
		random_inserted_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_real_body_states.WriteToFile(0.0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
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
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_complex);
	/** Initialize normal direction of the inserted body. */
	solid_dynamics::NormalDirectionSummation 	get_inserted_body_normal(inserted_body_complex);
	/** Corrected strong configuration for the elastic  insertbody. */
	solid_dynamics::CorrectConfiguration 		inserted_body_corrected_configuration_in_strong_form(inserted_body_inner);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_density(water_block_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationFirstHalf
		pressure_relaxation_first_half(water_block_complex);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAcceleration 	viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityFormulation 	transport_velocity_formulation(water_block_complex);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInFluidField 	compute_vorticity(water_block_inner);
	/** Inflow boundary condition. */
	ParabolicInflow		parabolic_inflow(water_block, new InflowBuffer(water_block, "Buffer"));
	/**
	 * @brief Algorithms of FSI.
	 */
	 /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidPressureForceOnSolid 	fluid_pressure_force_on_inserted_body(inserted_body_contact);
	solid_dynamics::FluidViscousForceOnSolid 	fluid_viscous_force_on_inserted_body(inserted_body_contact);
	/** Compute the average velocity of the insert body. */
	solid_dynamics::AverageVelocityAndAcceleration 		average_velocity_and_acceleration(inserted_body);
	/**
	 * @brief Algorithms of solid dynamics.
	 */
	 /** Compute time step size of elastic solid. */
	solid_dynamics::AcousticTimeStepSize 	inserted_body_computing_time_step_size(inserted_body);
	/** Stress relaxation for the inserted body. */
	solid_dynamics::StressRelaxationFirstHalf 	inserted_body_stress_relaxation_first_half(inserted_body_inner);
	solid_dynamics::StressRelaxationSecondHalf 	inserted_body_stress_relaxation_second_half(inserted_body_inner);
	/** Constrain region of the inserted body. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_beam_base(inserted_body, new BeamBase(inserted_body, "BeamBase"));
	/** Update norm .*/
	solid_dynamics::UpdateElasticNormalDirection 	inserted_body_update_normal(inserted_body);
	/**
	 * @brief Write observation data into files.
	 */
	WriteTotalViscousForceOnSolid 		write_total_viscous_force_on_inserted_body(in_output, inserted_body);

	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::pos_n_>
		write_beam_tip_displacement("Displacement", in_output, beam_observer_contact);
	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::vel_n_>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);

	/**
	 * @brief Pre-simulation.
	 */
	 /** Using relaxed particle distribution if needed. */
	if (system.reload_particles_) {
		unique_ptr<ReadReloadParticle>	
			reload_insert_body_particles(new ReadReloadParticle(in_output, { inserted_body }, { "InsertedBody" }));
		reload_insert_body_particles->ReadFromFile();
	}
	/** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	  * but before the configuration build up. */
	periodic_condition.parallel_exec();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
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
		inserted_body->updateCellLinkedList();
		water_block->updateCellLinkedList();
		periodic_condition.parallel_exec();
		/** one need update configuration after periodic condition. */
		water_block_complex->updateConfiguration();
		inserted_body_contact->updateConfiguration();
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
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time) {
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_viscous_force_on_inserted_body.parallel_exec();
			/** Update normal direction on elastic body.*/
			inserted_body_update_normal.parallel_exec();
			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				/** Fluid pressure relaxation, first half. */
				pressure_relaxation_first_half.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_inserted_body.parallel_exec();
				/** Fluid pressure relaxation, second half. */
				pressure_relaxation_second_half.parallel_exec(dt);

				/** Solid dynamics. */
				inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt) {

					dt_s = inserted_body_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;
					inserted_body_stress_relaxation_first_half.parallel_exec(dt_s);
					constrain_beam_base.parallel_exec();
					inserted_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
					inner_ite_dt_s++;
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
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

			water_block->updateCellLinkedList();
			periodic_condition.parallel_exec();
			water_block_complex->updateConfiguration();
			/** one need update configuration after periodic condition. */
			inserted_body->updateCellLinkedList();
			inserted_body_contact->updateConfiguration();
			/** write run-time observation into file */
			write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		write_total_viscous_force_on_inserted_body.WriteToFile(GlobalStaticVariables::physical_time_);
		fluid_observer_contact->updateConfiguration();
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

