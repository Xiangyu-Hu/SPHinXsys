/**
 * @file 	Dambreak.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;					/**< Tank length. */
Real DH = 5.366;					/**< Tank height. */
Real LL = 2.0;						/**< Liquid column length. */
Real LH = 1.0;						/**< Liquid column height. */
Real particle_spacing_ref = 0.025;	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;						 /**< Reference density of fluid. */
Real gravity_g = 1.0;					 /**< Gravity force of fluid. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin:
Vec2d water_block_translation = water_block_halfsize;
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transform2d(outer_wall_translation), outer_wall_halfsize);
		subtract<TransformShape<GeometricShapeBox>>(Transform2d(inner_wall_translation), inner_wall_halfsize);
	}
}; 
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	sph_system.handleCommandlineOptions(ac, av);
	/** I/O environment. */
	InOutput in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(
		sph_system, makeShared<TransformShape<GeometricShapeBox>>(
						Transform2d(water_block_translation), water_block_halfsize, "WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

	ObserverBody fluid_observer(sph_system, "FluidObserver");
	fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	Gravity gravity(Vecd(0.0, -gravity_g));
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	TimeStepInitialization fluid_step_initialization(water_block, gravity);
	fluid_dynamics::DensitySummationFreeSurfaceComplex fluid_density_by_summation(water_block_complex);
	fluid_dynamics::AdvectionTimeStepSize fluid_advection_time_step(water_block, U_max);
	fluid_dynamics::AcousticTimeStepSize fluid_acoustic_time_step(water_block);
	fluid_dynamics::PressureRelaxationRiemannWithWall fluid_pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall fluid_density_relaxation(water_block_complex);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	RestartIO restart_io(in_output, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<TotalMechanicalEnergy>>
		write_water_mechanical_energy(in_output, water_block, gravity);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_boundary_normal_direction.parallel_exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block.updateCellLinkedList();
		water_block_complex.updateConfiguration();
		fluid_observer_contact.updateConfiguration();
	}
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
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_fluid_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile();
	write_water_mechanical_energy.writeToFile(number_of_iterations);
	write_recorded_water_pressure.writeToFile(number_of_iterations);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** outer loop for dual-time criteria time-stepping. */
			time_instance = tick_count::now();
			fluid_step_initialization.parallel_exec();
			Real Dt = fluid_advection_time_step.parallel_exec();
			fluid_density_by_summation.parallel_exec();
			interval_computing_time_step += tick_count::now() - time_instance;

			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** inner loop for dual-time criteria time-stepping.  */
				fluid_pressure_relaxation.parallel_exec(dt);
				fluid_density_relaxation.parallel_exec(dt);
				dt = fluid_acoustic_time_step.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_fluid_pressure_relaxation += tick_count::now() - time_instance;

			/** screen output, write body reduced values and restart files  */
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.restart_step_)
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
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
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
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
			  << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
			  << interval_updating_configuration.seconds() << "\n";

	if (sph_system.restart_step_ == 0)
	{
		write_water_mechanical_energy.newResultTest();
		write_recorded_water_pressure.newResultTest();
	}

	return 0;
};
