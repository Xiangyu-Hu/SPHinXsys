/**
 * @file 	2d_flow_around_cylinder.cpp
 * @brief 	This is the benchmark test for the wall modeling of viscous flow.
 * @details We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu
 * @version 0.1
 */
#include "sphinxsys.h"

/** case file to setup the test case */
#include "2d_flow_around_cylinder.h"

using namespace SPH;

int main(int ac, char* av[])
{
	/** Build up -- a SPHSystem -- */
	SPHSystem system(Vec2d(- DL_sponge, - DH_sponge), Vec2d(DL, DH + DH_sponge), particle_spacing_ref);
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
	WaterBlock* water_block = new WaterBlock(system, "WaterBody", 0);
	WaterMaterial* water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Creating the cylinder.
	 */
	Cylinder* cylinder = new Cylinder(system, "Cylinder", 1);
	SolidParticles 	cylinder_particles(cylinder);
	/**
	 * @brief 	Particle and body creation of fluid observers.
	 */
	FluidObserver* fluid_observer = new FluidObserver(system, "FluidObserver", 0);
	BaseParticles	flow_observer_particles(fluid_observer);
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToVtu 				write_real_body_states(in_output, system.real_bodies_);
	WriteRestart						write_restart_files(in_output, system.real_bodies_);
	ReadRestart							read_restart_files(in_output, system.real_bodies_);
	/** body topology */
	SPHBodyComplexRelation* water_block_complex = new SPHBodyComplexRelation(water_block, {cylinder });
	SPHBodyContactRelation* cylinder_contact = new SPHBodyContactRelation(cylinder, { water_block });
	SPHBodyContactRelation* fluid_observer_contact = new SPHBodyContactRelation(fluid_observer, { water_block });

	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_) 
	{
		/** body topology only for particle realxation */
		SPHBodyInnerRelation* cylinder_inner = new SPHBodyInnerRelation(cylinder);
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		/** Random reset the insert body particle position. */
		RandomizePartilePosition  random_inserted_body_particles(cylinder);
		/** Write the body state to Vtu file. */
		WriteBodyStatesToVtu 		write_inserted_body_to_vtu(in_output, { cylinder });
		/** Write the particle reload files. */
		WriteReloadParticle 		write_particle_reload_files(in_output, { cylinder });

		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner);
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
		std::cout << "The physics relaxation process of the cylinder finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.WriteToFile(0.0);
		return 0;
	}
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList 	periodic_condition_x(water_block, 0);
	/** Periodic BCs in y direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList 	periodic_condition_y(water_block, 1);

	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_density(water_block_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationFirstHalf pressure_relaxation_first_half(water_block_complex);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann pressure_relaxation_second_half(water_block_complex);
	/** Computing viscous acceleration with wall model. */
	fluid_dynamics::ViscousAccelerationWallModel  viscous_acceleration_wall_modeling(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityFormulation 	transport_velocity_formulation(water_block_complex);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInFluidField 	compute_vorticity(water_block_complex->InnerRelation());
	/** freestream boundary condition. */
	FreeStreamCondition freestream_condition(water_block, new FreeStreamBuffer(water_block, "FreestreamBuffer"));
	/**
	 * @brief Algorithms of FSI.
	 */
	 /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidPressureForceOnSolid 	fluid_pressure_force_on_inserted_body(cylinder_contact);
	/** Computing viscous force acting on wall with wall model. */
	solid_dynamics::FluidViscousForceOnSolidWallModel
		            fluid_viscous_force_on_inserted_body_wall_modeling(cylinder_contact,
					&viscous_acceleration_wall_modeling);
	/**
	 * @brief Write observation data into files.
	 */
	WriteTotalViscousForceOnSolid 		write_total_viscous_force_on_inserted_body(in_output, cylinder);
	WriteTotalForceOnSolid              write_total_force_on_inserted_body(in_output, cylinder);

	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::vel_n_>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);

	/**
	 * @brief Pre-simulation.
	 */
	 /** Using relaxed particle distribution if needed. */
	if (system.reload_particles_) {
		unique_ptr<ReadReloadParticle>	
			reload_insert_body_particles(new ReadReloadParticle(in_output, { cylinder }, { "InsertedBody" }));
		reload_insert_body_particles->ReadFromFile();
	}
	/** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	  * but before the configuration build up. */
	periodic_condition_x.update_cell_linked_list_.parallel_exec();
	periodic_condition_y.update_cell_linked_list_.parallel_exec();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** initialize surface normal direction for the insert body. */
	cylinder_particles.initializeNormalDirectionFromGeometry();

	/**
	 * @brief The time stepping starts here.
	 */
	if (system.restart_step_ != 0) {
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(system.restart_step_);
		cylinder->updateCellLinkedList();
		water_block->updateCellLinkedList();
		periodic_condition_x.parallel_exec();
		periodic_condition_y.parallel_exec();
		/** one need update configuration after periodic condition. */
		water_block_complex->updateConfiguration();
		cylinder_contact->updateConfiguration();
	}
	/** first output*/
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 200.0;			    /**< End time. */
	Real D_Time = End_Time / 200.0;	/**< time stamps for output. */
	Real Dt = 0.0;					    /**< Default advection time step sizes for fluid. */
	Real dt = 0.0; 					    /**< Default acoustic time step sizes for fluid. */
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
		    viscous_acceleration_wall_modeling.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_viscous_force_on_inserted_body_wall_modeling.parallel_exec();
			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				/** Fluid pressure relaxation, first half. */
				pressure_relaxation_first_half.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_inserted_body.parallel_exec();
				/** Fluid pressure relaxation, second half. */
				pressure_relaxation_second_half.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				freestream_condition.parallel_exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					write_restart_files.WriteToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_condition_x.bounding_.parallel_exec();
			periodic_condition_y.bounding_.parallel_exec();
			water_block->updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.parallel_exec();
			periodic_condition_y.update_cell_linked_list_.parallel_exec();
			/** one need update configuration after periodic condition. */
			water_block_complex->updateConfiguration();
			cylinder_contact->updateConfiguration();

		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		write_total_viscous_force_on_inserted_body.WriteToFile(GlobalStaticVariables::physical_time_);
		write_total_force_on_inserted_body.WriteToFile(GlobalStaticVariables::physical_time_);
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

