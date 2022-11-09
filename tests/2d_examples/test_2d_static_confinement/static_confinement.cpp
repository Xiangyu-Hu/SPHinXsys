/**
 * @file 	static_confinement.cpp
 * @brief 	2D dambreak example in which the solid wall boundary are static confinement.
 * @details This is the one of the basic test cases.
 * @author 	Xiangyu Hu
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 5.366;			  /**< Tank length. */
Real DH = 5.366;			  /**< Tank height. */
Real LL = 2.0;				  /**< Liquid column length. */
Real LH = 1.0;				  /**< Liquid column height. */
Real resolution_ref = 0.025;  /**< Global reference resolution. */
Real BW = resolution_ref * 4; /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// Observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;						 /**< Reference density of fluid. */
Real gravity_g = 1.0;					 /**< Gravity force of fluid. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	// geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.0, 0.0));
	water_block_shape.push_back(Vecd(0.0, LH));
	water_block_shape.push_back(Vecd(LL, LH));
	water_block_shape.push_back(Vecd(LL, 0.0));
	water_block_shape.push_back(Vecd(0.0, 0.0));
	return water_block_shape;
}
/**
 * @brief create wall shape
 */
std::vector<Vecd> createWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

	return inner_wall_shape;
}

/** create a structure shape */
std::vector<Vecd> createStructureShape()
{
	// geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.5 * DL, 0.05 * DH));
	water_block_shape.push_back(Vecd(0.5 * DL + 0.5 * LL, 0.05 * DH + 0.5 * LH));
	water_block_shape.push_back(Vecd(0.5 * DL + 0.5 * LL, 0.05 * DH));
	water_block_shape.push_back(Vecd(0.5 * DL, 0.05 * DH));
	return water_block_shape;
}
/**
 *@brief 	Fluid body shape definition.
 */
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
	}
};
/**
 * @brief 	wall and structure surface definition.
 */
class WallAndStructure : public MultiPolygonShape
{
public:
	explicit WallAndStructure(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createStructureShape(), ShapeBooleanOps::sub);
	}
};
/**
 * @brief 	Main program starts here.
 */
int main()
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	/** output environment. */
	IOEnvironment io_environment(sph_system);
	/**
	 * @brief Material property, particles and body creation of fluid.
	 */
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	/**
	 * @brief 	Particle and body creation of fluid observer.
	 */
	ObserverBody fluid_observer(sph_system, "FluidObserver");
	fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	/** topology */
	InnerRelation water_block_inner(water_block);
	ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	/** Define external force. */
	Gravity gravity(Vecd(0.0, -gravity_g));
	/**
	 * @brief 	Methods used for time stepping.
	 */
	/** Initialize particle acceleration. */
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(water_block_inner);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using position verlet time stepping. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> pressure_relaxation(water_block_inner);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> density_relaxation(water_block_inner);
	/** Confinement condition for wall and structure. */
	NearShapeSurface near_surface(water_block, makeShared<WallAndStructure>("WallAndStructure"));
	fluid_dynamics::StaticConfinement confinement_condition(near_surface);
	update_density_by_summation.post_processes_.push_back(&confinement_condition.density_summation_);
	pressure_relaxation.post_processes_.push_back(&confinement_condition.pressure_relaxation_);
	density_relaxation.post_processes_.push_back(&confinement_condition.density_relaxation_);
	/**
	 * @brief Output.
	 */
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
	/** output the observed data from fluid body. */
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);

	/** Pre-simulation*/
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();

	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.writeToFile(0);
	write_recorded_water_pressure.writeToFile(0);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	Real end_time = 20.0;		/**< End time. */
	Real output_interval = 0.1; /**< Time stamps for output of body states. */
	Real dt = 0.0;				/**< Default acoustic time step sizes. */
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
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			interval_computing_time_step += tick_count::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
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
					write_recorded_water_pressure.writeToFile(number_of_iterations);
				}
				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();
			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_inner.updateConfiguration();
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
	write_recorded_water_pressure.newResultTest();

	return 0;
}
