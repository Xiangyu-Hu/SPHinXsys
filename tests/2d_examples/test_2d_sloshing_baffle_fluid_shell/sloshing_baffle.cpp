/**
* @file 	sloshing-baffle.cpp
* @brief 	Fluid-shell interaction in sloshing flow. 
* @details 	Here, the first fluid-shell interaction test is presented. 
* @author 	Chi Zhang
*/
#include "sphinxsys.h" 
#include "case.h"
using namespace SPH;

int main(int ac, char *av[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	IOEnvironment io_environment(sph_system);
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av);
#endif
	/** Set the starting time to zero. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** The water block, body, material and particles container. */
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	 /** The wall boundary, body and particles container. */
	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
	/** The baffle, body and particles container. */
	SolidBody shell_baffle(sph_system, makeShared<DefaultShape>("ShellBaffle"));
	shell_baffle.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	shell_baffle.generateParticles<ShellBaffleParticleGenerator>();
	shell_baffle.addBodyStateForRecording<Vecd>("PseudoNormal");
	/** @brief 	Particle and body creation of baffle observer.*/
	ObserverBody baffle_disp_observer(sph_system, "BaffleDispObserver");
	baffle_disp_observer.generateParticles<ObserverParticleGenerator>(baffle_disp_probe_location);
	/** Pressure probe on Flap. */
	ObserverBody baffle_pressure_observer(sph_system, "BafflePressureObserver");
	baffle_pressure_observer.generateParticles<ObserverParticleGenerator>(baffle_pressure_probe_location);
	/** @brief 	Particle and body creation of fluid observer.*/
	ObserverBody fluid_observer(sph_system, "Fluidobserver");
	fluid_observer.generateParticles<ObserverParticleGenerator>(fluid_pressure_probe_location);
	/** topology */
	InnerRelation 	water_block_inner(water_block);
	InnerRelation 	baffle_inner(shell_baffle);
	ContactRelation water_shell_contact(water_block, {&shell_baffle});
	ContactRelation water_wall_contact(water_block, {&wall_boundary});

	ComplexRelation water_wall_complex(water_block_inner, water_wall_contact);
	ComplexRelation water_shell_complex(water_block_inner, water_shell_contact);

	ContactRelation baffle_water_contact(shell_baffle, { &water_block });
	ContactRelation observer_contact_with_water(fluid_observer, { &water_block });
	ContactRelation observer_contact_with_baffle(baffle_disp_observer, { &shell_baffle});
	/** Density and wall norm. */
	SharedPtr<VariableGravity> gravity_ptr = makeShared<VariableGravity>();
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	/** Algorithm for Fluid dynamics. */
	SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block, gravity_ptr);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_fluid_density_by_summation(water_shell_contact, water_wall_complex);
	Dynamics1Level<fluid_dynamics::FluidShellandWallIntegration1stHalfRiemann> fluid_pressure_relaxation(water_wall_contact, water_shell_complex);
	Dynamics1Level<fluid_dynamics::FluidShellandWallIntegration2ndHalfRiemann> fluid_density_relaxation(water_wall_contact, water_shell_complex);
	/** Algorithms for solid. */
	ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell_baffle);
	InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(baffle_inner);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(baffle_inner);
	Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(baffle_inner);
	/** FSI */
	InteractionDynamics<solid_dynamics::FluidViscousForceOnShell> viscous_force_on_shell(baffle_water_contact);
	InteractionDynamics<solid_dynamics::FluidForceOnShellUpdate> fluid_force_on_shell_update(baffle_water_contact, viscous_force_on_shell);
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell_baffle);
	/** constraint and damping */
	BoundaryGeometry shell_boundary_geometry(shell_baffle, "BoundaryGeometry");
	SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion, BoundaryGeometry> baffle_constrain(shell_boundary_geometry);
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
		baffle_position_damping(0.2, baffle_inner, "Velocity", physical_viscosity);
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
		baffle_rotation_damping(0.2, baffle_inner, "AngularVelocity", physical_viscosity);
	/**
	 * @brief Output.
	 */
	BodyStatesRecordingToPlt write_real_body_states_to_plt(io_environment, sph_system.real_bodies_);
	BodyStatesRecordingToVtp write_real_body_states_to_vtp(io_environment, sph_system.real_bodies_);
	/** Output the observed displacement of baffle. */
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_baffle_displacement("Position", io_environment, observer_contact_with_baffle);
	/** Output the observed pressure of fluid. */
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_fluid_pressure_wall("Pressure", io_environment, observer_contact_with_water);
	/** WaveProbes. */
	/** #1.*/
	BodyRegionByCell wave_probe_buffer_no_1(water_block, makeShared<MultiPolygonShape>(CreateWaveProbeShape1(), "WaveProbe_01"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		wave_probe_1(io_environment, wave_probe_buffer_no_1);
	/** #2.*/
	BodyRegionByCell wave_probe_buffer_no_2(water_block, makeShared<MultiPolygonShape>(CreateWaveProbeShape2(), "WaveProbe_02"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		wave_probe_2(io_environment, wave_probe_buffer_no_2);
	/** #3.*/
	BodyRegionByCell wave_probe_buffer_no_3(water_block, makeShared<MultiPolygonShape>(CreateWaveProbeShape3(), "WaveProbe_03"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		wave_probe_3(io_environment, wave_probe_buffer_no_3);
	/**
	 * @brief The time stepping starts here.
	 */
	 /** Prepare quantities will be used once only and initial condition.*/
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_boundary_normal_direction.parallel_exec();
	shell_corrected_configuration.parallel_exec();
	/** Initial output. */
	//write_real_body_states_to_plt.writeToFile(0);
	write_real_body_states_to_vtp.writeToFile(0);
	wave_probe_1.writeToFile(0);
	wave_probe_2.writeToFile(0);
	wave_probe_3.writeToFile(0);
	write_baffle_displacement.writeToFile(0);
	write_fluid_pressure_wall.writeToFile(0);
	/** Time parameters. */
	int number_of_iterations = 0;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 5.0;			/**< End time. */
	Real D_Time = End_Time /100;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes. */
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
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
			/** Acceleration due to viscous force and gravity. */
			gravity_ptr->UpdateAcceleration();
			fluid_step_initialization.parallel_exec();
			update_fluid_density_by_summation.parallel_exec();

			viscous_force_on_shell.parallel_exec();

			Dt = 0.5 * fluid_advection_time_step.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** Fluid relaxation and force computaton. */
				fluid_pressure_relaxation.parallel_exec(dt);
				fluid_force_on_shell_update.parallel_exec();
				fluid_density_relaxation.parallel_exec(dt);

				/** Solid dynamics time stepping. */
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					dt_s = shell_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;
					shell_stress_relaxation_first.parallel_exec(dt_s);
					baffle_constrain.parallel_exec();
					shell_stress_relaxation_second.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				dt = fluid_acoustic_time_step.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedList();
			shell_baffle.updateCellLinkedList();

			water_block_inner.updateConfiguration();
			water_wall_contact.updateConfiguration();
			water_shell_contact.updateConfiguration();
			baffle_water_contact.updateConfiguration();
			observer_contact_with_water.updateConfiguration();

			wave_probe_1.writeToFile(number_of_iterations);
			wave_probe_2.writeToFile(number_of_iterations);
			wave_probe_3.writeToFile(number_of_iterations);
			write_baffle_displacement.writeToFile(number_of_iterations);
			write_fluid_pressure_wall.writeToFile(number_of_iterations);

		}
		write_real_body_states_to_vtp.writeToFile();

		tick_count t2 = tick_count::now();
	
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	if (sph_system.generate_regression_data_)
	{
		write_baffle_displacement.generateDataBase(1.0e-2);
	}
	else
	{
		write_baffle_displacement.newResultTest();
	}

	return 0;
}
