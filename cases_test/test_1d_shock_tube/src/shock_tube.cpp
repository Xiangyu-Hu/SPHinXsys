/**
 * @file 	shock_tube.cpp
 * @brief 	This is a test to show the sod shock tube.
 * @details We consider sod shock tube about two waves.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "shock_tube.h"
#include "sphinxsys.h"
using namespace SPH;
/**
 * @brief 	Main program starts here.
 */
int main(int ac, char* av[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	//handle command line arguments
	sph_system.handleCommandlineOptions(ac, av);
	/** output environment. */
	In_Output in_output(sph_system);
	/**
	 * @brief Material property, partilces and body creation of water.
	 */
	WaveBlock* wave_block = new WaveBlock(sph_system, "WaveBody");
	WaveMaterial* wave_material = new WaveMaterial();
	CompressibleFluidParticles 	wave_particles(wave_block, wave_material);
	wave_particles.addAVariableToWrite<indexScalar, Real>("TotalEnergy");
	/** topology */
	BaseBodyRelationInner* wave_block_inner = new BodyRelationInner(wave_block);
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Initial momentum, pressure and energy field */
	WavesInitialCondition setup_waves_initial_parameters(wave_block);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	eulerian_fluid_dynamics::CompressibleFlowTimeStepInitialization		initialize_wave_step(wave_block);
	/** Periodic BCs in y direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList 	periodic_condition_y(wave_block, yAxis);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Time step size with considering sound wave speed. */
	eulerian_fluid_dynamics::AcousticTimeStepSize get_wave_time_step_size(wave_block);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	/** Here, we can use HLLC Riemann solver for pressure relaxation and density and energy relaxation  */
	eulerian_fluid_dynamics::PressureRelaxationHLLCRiemannInner pressure_relaxation(wave_block_inner);
	eulerian_fluid_dynamics::DensityAndEnergyRelaxationHLLCRiemannInner density_and_energy_relaxation(wave_block_inner);
	/**
	 * @brief Output.
	 */
	 /** Output the body states. */
	BodyStatesRecordingToPlt 		body_states_recording(in_output, sph_system.real_bodies_);
	/** Write the particle reload files. */
	ReloadParticleIO 		write_particle_reload_files(in_output, { wave_block });
	/** Output the body states for restart simulation. */
	RestartIO		restart_io(in_output, sph_system.real_bodies_);
	/**
	 * @brief Setup geomtry and initial conditions
	 */
	setup_waves_initial_parameters.exec();
	sph_system.initializeSystemCellLinkedLists();
	periodic_condition_y.update_cell_linked_list_.parallel_exec();
	sph_system.initializeSystemConfigurations();
	/**
	 * @brief The time stepping starts here.
	 */
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		wave_block->updateCellLinkedList();
		periodic_condition_y.update_cell_linked_list_.parallel_exec();
		wave_block_inner->updateConfiguration();
	}
	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 0.2; 	/**< End time. */
	Real D_Time = 0.01;		/**< Time stamps for output of body states. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * @brief 	Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_wave_step.parallel_exec();
			dt = get_wave_time_step_size.parallel_exec();
			/** Dynamics including pressure and density and energy relaxation. */
			integration_time += dt;
			pressure_relaxation.parallel_exec(dt);
			density_and_energy_relaxation.parallel_exec(dt);
			GlobalStaticVariables::physical_time_ += dt;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0) {
					restart_io.writeToFile(number_of_iterations);
				}
			}
			number_of_iterations++;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << endl;

	write_particle_reload_files.writeToFile();
	return 0;
}
