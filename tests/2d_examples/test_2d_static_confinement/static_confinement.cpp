/**
 * @file 	static_confinement.cpp
 * @brief 	2D dambreak example in which the solid wall boundary are static confinement.
 * @details This is the one of the basic test cases.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;			  /**< Tank length. */
Real DH = 5.366;			  /**< Tank height. */
Real LL = 2.0;				  /**< Liquid column length. */
Real LH = 1.0;				  /**< Liquid column height. */
Real resolution_ref = 0.025;  /**< Global reference resolution. */
Real BW = resolution_ref * 4; /**< Extending width for BCs. */
// Observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;						 /**< Reference density of fluid. */
Real gravity_g = 1.0;					 /**< Gravity force of fluid. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
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
/** create wall shape */
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
//----------------------------------------------------------------------
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Shape for wall boundary and a structure.
//----------------------------------------------------------------------
class WallAndStructure : public MultiPolygonShape
{
public:
	explicit WallAndStructure(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createStructureShape(), ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	ObserverBody fluid_observer(sph_system, "FluidObserver");
	fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_block);
	ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** Initialize particle acceleration. */
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, gravity_ptr);
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
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(io_environment, water_block, gravity_ptr);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	Real end_time = 20.0;		/**< End time. */
	Real output_interval = 0.1; /**< Time stamps for output of body states. */
	Real dt = 0.0;				/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	TimeInterval interval_computing_time_step;
	TimeInterval interval_computing_pressure_relaxation;
	TimeInterval interval_updating_configuration;
	TickCount time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile(0);
	write_water_mechanical_energy.writeToFile(0);
	write_recorded_water_pressure.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = TickCount::now();
			initialize_a_fluid_step.exec();
			Real Dt = get_fluid_advection_time_step_size.exec();
			update_density_by_summation.exec();
			interval_computing_time_step += TickCount::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation.exec(dt);
				density_relaxation.exec(dt);
				dt = get_fluid_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += TickCount::now() - time_instance;

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
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = TickCount::now();
			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_inner.updateConfiguration();
			fluid_observer_contact.updateConfiguration();
			interval_updating_configuration += TickCount::now() - time_instance;
		}

		TickCount t2 = TickCount::now();
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
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
