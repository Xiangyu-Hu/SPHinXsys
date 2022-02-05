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

int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	/** I/O environment. */
	In_Output in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	WaterBlock water_block(sph_system, "WaterBody");
	FluidParticles water_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f));

	AirBlock air_block(sph_system, "AirBody");
	FluidParticles air_particles(air_block, makeShared<WeaklyCompressibleFluid>(rho0_a, c_f));

	WallBoundary wall_boundary(sph_system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	wall_particles.addAVariableToWrite<Vecd>("NormalDirection");

	ObserverBody fluid_observer(sph_system, "Fluidobserver");
	ObserverParticles observer_particles(fluid_observer, makeShared<ObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexBodyRelation water_air_complex(water_block, {&air_block});
	BodyRelationContact water_wall_contact(water_block, {&wall_boundary});
	ComplexBodyRelation air_water_complex(air_block, {&water_block});
	BodyRelationContact air_wall_contact(air_block, {&wall_boundary});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block, &air_block});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Define external force. */
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_water_step(water_block, gravity);
	TimeStepInitialization initialize_a_air_step(air_block, gravity);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex
		update_water_density_by_summation(water_air_complex.inner_relation_, water_wall_contact);
	fluid_dynamics::DensitySummationComplex
		update_air_density_by_summation(air_water_complex, air_wall_contact);
	fluid_dynamics::TransportVelocityCorrectionComplex
		air_transport_correction(air_water_complex, air_wall_contact);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_water_advection_time_step_size(water_block, U_max);
	fluid_dynamics::AdvectionTimeStepSize get_air_advection_time_step_size(air_block, U_max);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_water_time_step_size(water_block);
	fluid_dynamics::AcousticTimeStepSize get_air_time_step_size(air_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	fluid_dynamics::PressureRelaxationRiemannWithWall
		water_pressure_relaxation(water_air_complex.inner_relation_, water_wall_contact);
	fluid_dynamics::DensityRelaxationRiemannWithWall
		water_density_relaxation(water_air_complex.inner_relation_, water_wall_contact);
	/** Extend Pressure relaxation is used for air. */
	fluid_dynamics::ExtendMultiPhasePressureRelaxationRiemannWithWall
		air_pressure_relaxation(air_water_complex, air_wall_contact, 2.0);
	fluid_dynamics::MultiPhaseDensityRelaxationRiemannWithWall
		air_density_relaxation(air_water_complex, air_wall_contact);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations 
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(in_output, sph_system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<TotalMechanicalEnergy>>
		write_water_mechanical_energy(in_output, water_block, gravity);
	/** output the observed data from fluid body. */
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_pressure("Pressure", in_output, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromShapeAndOp("InnerWall");
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		water_wall_contact.updateConfiguration();
		air_water_complex.updateConfiguration();
		air_wall_contact.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
 	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.writeToFile(0);
	write_recorded_pressure.writeToFile(0);
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 20.0; /**< End time. */
	Real D_Time = 0.1;	  /**< Time stamps for output of body states. */
	Real dt = 0.0;		  /**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			initialize_a_water_step.parallel_exec();
			initialize_a_air_step.parallel_exec();

			Real Dt_f = get_water_advection_time_step_size.parallel_exec();
			Real Dt_a = get_air_advection_time_step_size.parallel_exec();
			Real Dt = SMIN(Dt_f, Dt_a);

			update_water_density_by_summation.parallel_exec();
			update_air_density_by_summation.parallel_exec();

			air_transport_correction.parallel_exec(Dt);

			interval_computing_time_step += tick_count::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.parallel_exec();
				Real dt_a = get_air_time_step_size.parallel_exec();
				dt = SMIN(SMIN(dt_f, dt_a), Dt);

				water_pressure_relaxation.parallel_exec(dt);
				air_pressure_relaxation.parallel_exec(dt);

				water_density_relaxation.parallel_exec(dt);
				air_density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
				{
					write_water_mechanical_energy.writeToFile(number_of_iterations);
					write_recorded_pressure.writeToFile(number_of_iterations);
				}
				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();

			water_block.updateCellLinkedList();
			water_air_complex.updateConfiguration();
			water_wall_contact.updateConfiguration();

			air_block.updateCellLinkedList();
			air_water_complex.updateConfiguration();
			air_wall_contact.updateConfiguration();

			fluid_observer_contact.updateConfiguration();
			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}

	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
			  << interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
			  << interval_computing_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
			  << interval_updating_configuration.seconds() << "\n";

	write_water_mechanical_energy.newResultTest();
	write_recorded_pressure.newResultTest();

	return 0;
}
