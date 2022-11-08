/**
 * @file 	droplet.cpp
 * @brief 	A square droplet deforms to circle due to surface tension.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "droplet.h"
#include "sphinxsys.h"
using namespace SPH;
/**
 * @brief 	Main program starts here.
 */
int main()
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/**
	 * @brief Material property, particles and body creation of water.
	 */
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	water_block.addBodyStateForRecording<int>("SurfaceIndicator");
	/**
	 * @brief Material property, particles and body creation of air.
	 */
	FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_a, c_f, mu_a);
	air_block.generateParticles<ParticleGeneratorLattice>();
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	/** body topology */
	ComplexRelation water_air_complex(water_block, {&air_block});
	ContactRelation water_wall_contact(water_block, {&wall_boundary});
	ComplexRelation air_water_complex(air_block, {&water_block});
	ContactRelation air_wall_contact(air_block, {&wall_boundary});
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	/** Initialize particle acceleration. */
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block);
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex>
		update_water_density_by_summation(water_wall_contact, water_air_complex.inner_relation_);
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
		update_air_density_by_summation(air_wall_contact, air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex>
		air_transport_correction(air_wall_contact, air_water_complex);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step_size(air_block, U_max);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_time_step_size(air_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall>
		water_pressure_relaxation(water_wall_contact, water_air_complex.inner_relation_);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall>
		water_density_relaxation(water_wall_contact, water_air_complex.inner_relation_);
	/** Extend Pressure relaxation is used for air. */
	Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
		air_pressure_relaxation(air_wall_contact, air_water_complex, 2.0);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		air_density_relaxation(air_wall_contact, air_water_complex);
	/** Viscous acceleration. */
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> air_viscous_acceleration(air_wall_contact, air_water_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> water_viscous_acceleration(water_wall_contact, water_air_complex);
	/** Surface tension. */
	InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationInner> surface_detection(water_air_complex.inner_relation_);
	InteractionDynamics<fluid_dynamics::ColorFunctionGradientInner> color_gradient(water_air_complex.inner_relation_);
	InteractionDynamics<fluid_dynamics::ColorFunctionGradientInterpolationInner> color_gradient_interpolation(water_air_complex.inner_relation_);
	InteractionDynamics<fluid_dynamics::SurfaceTensionAccelerationInner> surface_tension_acceleration(water_air_complex.inner_relation_, 1.0);
	/**
	 * @brief Output.
	 */
	IOEnvironment io_environment(sph_system);
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	/** Pre-simulation*/
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_boundary_normal_direction.parallel_exec();
	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	Real end_time = 1.0;
	Real output_interval = 0.02;	 /**< Time stamps for output of body states. */
	Real dt = 0.0;		 /**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;

	/**
	 * @brief 	Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
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
			air_transport_correction.parallel_exec();

			air_viscous_acceleration.parallel_exec();
			water_viscous_acceleration.parallel_exec();

			surface_detection.parallel_exec();
			color_gradient.parallel_exec();
			color_gradient_interpolation.parallel_exec();
			surface_tension_acceleration.parallel_exec();

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

			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();

			water_block.updateCellLinkedListWithParticleSort(100);
			water_air_complex.updateConfiguration();
			water_wall_contact.updateConfiguration();

			air_block.updateCellLinkedListWithParticleSort(100);
			air_water_complex.updateConfiguration();
			air_wall_contact.updateConfiguration();

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

	return 0;
}
