/**
 * @file 	waterentry.cpp
 * @brief 	This is a cylinder fall into static water.
 */
#include "sphinxsys.h"
#include "tank_case.h"
 /*@brief Namespace cite here.
 */
using namespace SPH;



/*
Main program starts here.
*/
int main(int ac, char* av[])
{
	/* Build up -- a SPHSystem -- */
	SPHSystem system(system_domain_bounds, resolution_ref);
	// Tag for run particle relaxation for the initial body fitted distribution.
	system.setRunParticleRelaxation(false);
	// Tag for computation start with relaxed body fitted particles distribution.
	system.setReloadParticles(true);
	/* Tag for computation from restart files. 0: start with initial condition. */
	system.setRestartStep(0);
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);
	/* Output environment. */
	IOEnvironment in_output(system);


	/*
	@Brief creating body, materials and particles for the cylinder.
	*/
	SolidBody tank(system, makeShared<Tank>("Tank"));
	tank.defineParticlesAndMaterial<SolidParticles, Solid>();
	/*tank.defineBodyLevelSetShape()->writeLevelSet(in_output);*/
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? tank.generateParticles<ParticleGeneratorReload>(in_output, tank.getName())
		: tank.generateParticles<ParticleGeneratorLattice>();


	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<TwoPhaseDiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoWaterBodyMaterial>();
	//water_block.defineBodyLevelSetShape()->writeLevelSet(in_output);
	water_block.generateParticles<ParticleGeneratorLattice>();


	FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<TwoPhaseDiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoAirBodyMaterial>();
	//air_block.defineBodyLevelSetShape()->writeLevelSet(in_output);
	air_block.generateParticles<ParticleGeneratorLattice>();

	/*ObserverBody liquid_temperature_observer(system, "LiquidTemperatureObserver");
	liquid_temperature_observer.generateParticles<LiquidTemperatureObserverParticleGenerator>();

	ObserverBody gas_temperature_observer(system, "GasTemperatureObserver");
	gas_temperature_observer.generateParticles<GasTemperatureObserverParticleGenerator>();
	*/
	InnerRelation tank_inner(tank);
	//InnerRelation water_inner(water_block);
	//InnerRelation air_inner(air_block);

	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
		//SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
		//SimpleDynamics<RandomizeParticlePosition> random_air_particles(air_block);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_tank_to_vtp(in_output, { &tank });
		//BodyStatesRecordingToVtp write_water_to_vtp(in_output, { &water_block });
		//BodyStatesRecordingToVtp write_air_to_vtp(in_output, { &air_block });
		/** Write the particle reload files. */
		ReloadParticleIO write_tank_particle_reload_files(in_output, tank, "Tank");
		//ReloadParticleIO write_water_particle_reload_files(in_output, water_block, "WaterBody");
	//	ReloadParticleIO write_air_particle_reload_files(in_output, air_block, "AirBody");
		/** A  Physics relaxation step. */
		//relax_dynamics::RelaxationStepInner water_relaxation_step_inner(water_inner);
		//relax_dynamics::RelaxationStepInner air_relaxation_step_inner(air_inner);
		relax_dynamics::RelaxationStepInner tank_relaxation_step_inner(tank_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_tank_particles.parallel_exec(0.25);
		//random_water_particles.parallel_exec(0.25);
		//random_air_particles.parallel_exec(0.25);
		tank_relaxation_step_inner.SurfaceBounding().parallel_exec();
		//water_relaxation_step_inner.SurfaceBounding().parallel_exec();
		//air_relaxation_step_inner.SurfaceBounding().parallel_exec();
		write_tank_to_vtp.writeToFile(0);
		//write_water_to_vtp.writeToFile(0);
		//write_air_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			tank_relaxation_step_inner.parallel_exec();
			//water_relaxation_step_inner.parallel_exec();
			//air_relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
				//write_water_to_vtp.writeToFile(ite_p);
				//write_air_to_vtp.writeToFile(ite_p);
				write_tank_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of finish !" << std::endl;
		/** Output results. */
		write_tank_particle_reload_files.writeToFile(0);
		//	write_water_particle_reload_files.writeToFile(0);
		//	write_air_particle_reload_files.writeToFile(0);
		return 0;
	}


	ContactRelation water_block_contact(water_block, { &tank });
	ContactRelation air_block_contact(air_block, { &tank });
	//ContactRelation liquid_temperature_observer_contact(liquid_temperature_observer, { &water_block });
	//ContactRelation gas_temperature_observer_contact(gas_temperature_observer, { &air_block });
	ComplexRelation water_air_complex(water_block, { &air_block });
	ComplexRelation air_water_complex(air_block, { &water_block });


	/*
	@Brief define simple data file input and outputs functions.
	*/
	BodyStatesRecordingToVtp 			write_real_body_states(in_output, system.real_bodies_);
	RestartIO							restart_io(in_output, system.real_bodies_);
	//ObservedQuantityRecording<Real> write_temperature_liquid("Phi", in_output, liquid_temperature_observer_contact);
	//ObservedQuantityRecording<Real> write_temperature_gas("Phi", in_output, gas_temperature_observer_contact);

	//InteractionDynamics<solid_dynamics::CorrectConfiguration> 		tank_corrected_configuration(tank_complex.inner_relation_);
	SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(tank, "InnerWall");
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, makeShared<VariableGravity>());
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, makeShared<VariableGravity>());
	/* Fluid dynamics */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> fluid_density_by_summation(water_block_contact, water_air_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> air_density_by_summation(air_block_contact, air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(air_block_contact, air_water_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_water(water_block_contact, water_air_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_air(air_block_contact, air_water_complex);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> air_advection_time_step(air_block, U_g);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> air_acoustic_time_step(air_block);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> fluid_pressure_relaxation(water_block_contact, water_air_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> fluid_density_relaxation(water_block_contact, water_air_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
		air_pressure_relaxation(air_block_contact, air_water_complex, 2.0);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		air_density_relaxation(air_block_contact, air_water_complex);

	SimpleDynamics<ThermoAirBodyInitialCondition> thermo_air_initial_condition(air_block);
	SimpleDynamics<ThermoWaterBodyInitialCondition> thermo_water_initial_condition(water_block);

	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_water(water_block);
	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_air(air_block);

	ThermalRelaxationComplex thermal_relaxation_complex_water(water_air_complex);
	ThermalRelaxationComplex thermal_relaxation_complex_air(air_water_complex);

	BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("PorbeS1"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_1(in_output, probe_s1);
	BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_2(in_output, probe_s2);
	BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("PorbeS3"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_3(in_output, probe_s3);

	ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		water_average_temperature(in_output, water_block, "Phi");
	ReducedQuantityRecording<ReduceAverage<DiffusionReactionSpeciesSummation<FluidParticles, WeaklyCompressibleFluid>>>
		air_average_temperature(in_output, air_block, "Phi");
	ReducedQuantityRecording<ReduceDynamics<QuantityMoment<Real>>>
		air_rate_of_heat_transfer(in_output, air_block, "HeatFlux");
	ReducedQuantityRecording<ReduceDynamics<QuantityMoment<Real>>>
		water_rate_of_heat_transfer(in_output, water_block, "HeatFlux");
	ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_air_total_mass(in_output, air_block, "MassiveMeasure");
	ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Real>>> compute_water_total_mass(in_output, water_block, "MassiveMeasure");
	/**
 * @brief Pre-simulation.
 */
 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the tank. */
	inner_normal_direction.parallel_exec();
	/** computing linear reproducing configuration for the tank. */
	//tank_corrected_configuration.parallel_exec();
	write_real_body_states.writeToFile(0);
	thermo_water_initial_condition.parallel_exec();
	thermo_air_initial_condition.parallel_exec();
	probe_1.writeToFile(0);
	probe_2.writeToFile(0);
	probe_3.writeToFile(0);
	water_average_temperature.writeToFile(0);
	air_average_temperature.writeToFile(0);
	air_rate_of_heat_transfer.writeToFile(0);
	water_rate_of_heat_transfer.writeToFile(0);
	compute_water_total_mass.writeToFile(0);
	compute_air_total_mass.writeToFile(0);
	if (system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		water_block_contact.updateConfiguration();
		air_water_complex.updateConfiguration();
		air_block_contact.updateConfiguration();
	}

	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 6.0;			/**< End time. */
	Real D_Time = 0.1;	/**< time stamps for output. */
	//Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes for fluid. */

	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** outer loop for dual-time criteria time-stepping. */
			initialize_a_water_step.parallel_exec();
			initialize_a_air_step.parallel_exec();
			Real Dt_f = fluid_advection_time_step.parallel_exec();
			Real Dt_a = air_advection_time_step.parallel_exec();
			Real Dt = SMIN(Dt_f, Dt_a);
			fluid_density_by_summation.parallel_exec();
			air_density_by_summation.parallel_exec();
			air_transport_correction.parallel_exec();
			viscous_acceleration_air.parallel_exec();
			viscous_acceleration_water.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = fluid_acoustic_time_step.parallel_exec();
				Real dt_a = air_acoustic_time_step.parallel_exec();
				Real dt_thermal_water = get_thermal_time_step_water.parallel_exec();
				Real dt_thermal_air = get_thermal_time_step_air.parallel_exec();
				dt = SMIN(SMIN(dt_f, dt_thermal_water), SMIN(dt_thermal_air, dt_a), Dt);
				/* Fluid pressure relaxation */
				fluid_pressure_relaxation.parallel_exec(dt);
				air_pressure_relaxation.parallel_exec(dt);
				/* Fluid density relaxation */
				fluid_density_relaxation.parallel_exec(dt);
				air_density_relaxation.parallel_exec(dt);
				thermal_relaxation_complex_air.parallel_exec(dt);
				thermal_relaxation_complex_water.parallel_exec(dt);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

			}

			/** screen output, write body reduced values and restart files  */
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}

			number_of_iterations++;

			/** Update cell linked list and configuration. */

			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_contact.updateConfiguration();
			water_air_complex.updateConfiguration();

			air_block.updateCellLinkedListWithParticleSort(100);
			air_block_contact.updateConfiguration();
			air_water_complex.updateConfiguration();
			air_rate_of_heat_transfer.writeToFile();
			water_rate_of_heat_transfer.writeToFile();
			//liquid_temperature_observer_contact.updateConfiguration();
			//gas_temperature_observer_contact.updateConfiguration();
			int num = air_water_complex.getInnerRelation().base_particles_.pos_.size();
			int count_num = 0;
			for (size_t k = 0; k < air_water_complex.contact_configuration_.size(); ++k)
			{
				for (int i = 0; i != num; ++i)
				{
					Neighborhood& neighborhood = air_water_complex.contact_configuration_[k][i];
					if (neighborhood.current_size_ != 0)
					{
						count_num += 1;
					}
				}
			}
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::string output_folder_ = "./output";
				std::string filefullpath = output_folder_ + "/" + "number_of_contact_particles" + ".dat";
				std::ofstream out_file_(filefullpath.c_str(), std::ios::app);
				out_file_ << GlobalStaticVariables::physical_time_ << "  " << count_num << "  " << std::endl;
			}
		}
		tick_count t2 = tick_count::now();
		/** write run-time observation into file */

		write_real_body_states.writeToFile();
		probe_1.writeToFile();
		probe_2.writeToFile();
		probe_3.writeToFile();
		water_average_temperature.writeToFile();
		air_average_temperature.writeToFile();
		//write_temperature_liquid.writeToFile();
		//write_temperature_gas.writeToFile();



		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}
