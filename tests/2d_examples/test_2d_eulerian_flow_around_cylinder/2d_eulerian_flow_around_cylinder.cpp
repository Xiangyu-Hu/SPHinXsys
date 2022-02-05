/**
 * @file 	2d_eulerian_flow_around_cylinder.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder.
 * @details We consider a flow passing by a cylinder in 2D.
 * @author 	Zhentong Wang
 */
#include "sphinxsys.h"

 /** case file to setup the test case */
#include "2d_eulerian_flow_around_cylinder.h"

using namespace SPH;

int main(int ac, char* av[])
{
	/** Build up -- a SPHSystem -- */
	SPHSystem system(system_domain_bounds, resolution_ref);
	In_Output	in_output(system);
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
	/**
	 * @brief Creating body, materials and particles for a water block.
	 */
	WaterBlock water_block(system, "WaterBody");
	SharedPtr<ParticleGenerator> water_block__particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		water_block__particle_generator = makeShared<ParticleGeneratorReload>(in_output, water_block.getBodyName());
	WeaklyCompressibleFluidParticles 	fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f),water_block__particle_generator);
	fluid_particles.addAVariableToWrite<Real>("Pressure");
	fluid_particles.addAVariableToWrite<int>("SurfaceIndicator");
	/**
	 * @brief 	Creating the cylinder.
	 */
	Cylinder cylinder(system, "Cylinder");
	SharedPtr<ParticleGenerator> cylinder_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		cylinder_particle_generator = makeShared<ParticleGeneratorReload>(in_output, cylinder.getBodyName());
	SolidParticles 	cylinder_particles(cylinder, cylinder_particle_generator);

	/**
	 * @brief define simple data file input and outputs functions.
	 */
	BodyStatesRecordingToVtp 				write_real_body_states(in_output, system.real_bodies_);
	RestartIO							restart_io(in_output, system.real_bodies_);
	/** body topology */
	BodyRelationInner water_block_inner(water_block);
	BodyRelationContact water_cylinder_contact(water_block, { &cylinder });
	ComplexBodyRelation water_block_complex(water_block, { &cylinder });
	BodyRelationContact cylinder_contact (cylinder, { &water_block });

	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_)
	{
		//body topology for particle realxation 
		BodyRelationInner cylinder_inner(cylinder);
		/*
		 * @brief 	Methods used for particle relaxation.
		 */
		 /** Random reset the insert body particle position. */
		RandomizePartilePosition  random_inserted_body_particles(cylinder);
		RandomizePartilePosition  random_water_body_particles(water_block);
		/** Write the body state to Vtu file. */
		BodyStatesRecordingToVtp 		write_inserted_body_to_vtu(in_output, { &cylinder });
		BodyStatesRecordingToVtp 		write_water_body_to_vtu(in_output, { &water_block });
		/** Write the particle reload files. */
		ReloadParticleIO 		write_inserted_particle_reload_files(in_output, { &cylinder });
		ReloadParticleIO 		write_water_particle_reload_files(in_output, { &water_block });

		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(cylinder_inner, true);
		relax_dynamics::RelaxationStepComplex relaxation_step_complex(water_block_complex, true);
		/**
		  * @brief 	Particle relaxation starts here.
		  */
		random_inserted_body_particles.parallel_exec(0.25);
		random_water_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		relaxation_step_complex.surface_bounding_.parallel_exec();
		write_real_body_states.writeToFile(0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			relaxation_step_complex.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps for the inserted and water body N = " << ite_p << "\n";
				write_inserted_body_to_vtu.writeToFile(ite_p);
				write_water_body_to_vtu.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of the cylinder and water block finish !" << std::endl;

		/** Output results. */
		write_inserted_particle_reload_files.writeToFile(0);
		write_water_particle_reload_files.writeToFile(0);
		return 0;
	}

	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	eulerian_weakly_compressible_fluid_dynamics::EulerianFlowTimeStepInitialization 	initialize_a_fluid_step(water_block);
	/** Time step size with considering sound wave speed. */
	eulerian_weakly_compressible_fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	eulerian_weakly_compressible_fluid_dynamics::PressureRelaxationHLLCRiemannWithLimiterWithWall pressure_relaxation(water_block_complex);
	eulerian_weakly_compressible_fluid_dynamics::DensityAndEnergyRelaxationHLLCRiemannWithLimiterWithWall density_relaxation(water_block_complex);
	eulerian_weakly_compressible_fluid_dynamics::FreeSurfaceIndicationComplex surface_indicator(water_block_complex.inner_relation_, water_block_complex.contact_relation_);
	/** Computing viscous acceleration with wall model. */
	eulerian_weakly_compressible_fluid_dynamics::ViscousAccelerationWithWall  viscous_acceleration(water_block_complex);
	/** non_reflective boundary condition. */
	FarFieldBoundary variable_reset_in_boundary_condition(water_block_complex.inner_relation_);
	/**
	 * @brief Algorithms of FSI.
	 */
	 /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidForceOnSolidUpdateRiemannWithLimiterInEuler fluid_force_on_solid_update(cylinder_contact);
	/**
	 * @brief Write solid data into files.
	 */
	RegressionTestTimeAveraged<BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid>>
		write_total_viscous_force_on_inserted_body(in_output, cylinder);
	BodyReducedQuantityRecording<solid_dynamics::TotalForceOnSolid>
		write_total_force_on_inserted_body(in_output, cylinder);
	/** Output the maximum speed of the fluid body. */
	BodyReducedQuantityRecording<MaximumSpeed> write_maximum_speed(in_output, water_block);
	/**
	 * @brief Pre-simulation.
	 */
	 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** initialize surface normal direction for the insert body. */
	cylinder_particles.initializeNormalDirectionFromBodyShape();
	surface_indicator.parallel_exec();
	variable_reset_in_boundary_condition.parallel_exec();

	/**
	 * @brief The time stepping starts here.
	 */
	if (system.restart_step_ != 0) {
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		cylinder.updateCellLinkedList();
		water_block.updateCellLinkedList();
		water_block_complex.updateConfiguration();
		cylinder_contact.updateConfiguration();
	}
	/** first output*/
	write_real_body_states.writeToFile(0);

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 200.0;			    /**< End time. */
	Real D_Time = 1.0;	/**< time stamps for output. */
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
			dt = get_fluid_time_step_size.parallel_exec();
			viscous_acceleration.parallel_exec();
			/** FSI for viscous force. */
			fluid_force_on_solid_update.viscous_force_.parallel_exec();
			/** Fluid pressure relaxation, first half. */
			pressure_relaxation.parallel_exec(dt);
			/** FSI for pressure force. */
			fluid_force_on_solid_update.parallel_exec();
			/** Fluid pressure relaxation, second half. */
			density_relaxation.parallel_exec(dt);

			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
			variable_reset_in_boundary_condition.parallel_exec();

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(Real(number_of_iterations));
			}
			number_of_iterations++;
		}

		tick_count t2 = tick_count::now();
		write_real_body_states.writeToFile();
		write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
		write_total_force_on_inserted_body.writeToFile(number_of_iterations);
		write_maximum_speed.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	write_total_viscous_force_on_inserted_body.newResultTest();

	return 0;
}

