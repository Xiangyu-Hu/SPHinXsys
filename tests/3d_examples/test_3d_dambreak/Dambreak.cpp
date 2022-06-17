/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/

#include "sphinxsys.h" // SPHinXsys Library.

using namespace SPH;

//for geometry
Real resolution_ref = 0.05;	  //particle spacing
Real BW = resolution_ref * 4; //boundary width
Real DL = 5.366;			  //tank length
Real DH = 2.0;				  //tank height
Real DW = 0.5;				  //tank width
Real LL = 2.0;				  //liquid length
Real LH = 1.0;				  //liquid height
Real LW = 0.5;				  //liquid width
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));

//for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;

//	define the fluid body shape
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
		SimTK::Transform translation_water(halfsize_water);
		add<TransformShape<GeometricShapeBox>>(translation_water, halfsize_water);
	}
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
		Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
		SimTK::Transform  translation_wall(halfsize_inner);
		add<TransformShape<GeometricShapeBox>>(translation_wall, halfsize_outer);
		subtract<TransformShape<GeometricShapeBox>>(translation_wall, halfsize_inner);
	}
};

//	define an observer particle generator
class WaterObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit WaterObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		//add observation points
		positions_.push_back(Vecd(DL, 0.01, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.1, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.2, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.24, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.252, 0.5 * DW));
		positions_.push_back(Vecd(DL, 0.266, 0.5 * DW));
	}
};

//the main program
int main()
{

	//build up context -- a SPHSystem
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;

	//the water block
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles,WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	//the wall boundary
	SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

	ObserverBody fluid_observer(system, "Fluidobserver");
	fluid_observer.generateParticles<WaterObserverParticleGenerator>();

	/** topology */
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------
	//-------- common particle dynamics ----------------------------------------
	Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	TimeStepInitialization initialize_a_fluid_step(water_block, gravity);
	//-------- fluid dynamics --------------------------------------------------
	//evaluation of density by summation approach
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex);
	//time step size without considering sound wave speed
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);

	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	InOutput in_output(system);
	BodyStatesRecordingToVtp write_water_block_states(in_output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(in_output, system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	RegressionTestEnsembleAveraged<BodyReducedQuantityRecording<TotalMechanicalEnergy>>
		write_water_mechanical_energy(in_output, water_block, gravity);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);
	//-------------------------------------------------------------------
	//from here the time stepping begins
	//-------------------------------------------------------------------

	/**
	 * @brief Setup geometrics and initial conditions
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_normal_direction.parallel_exec();
	/**
	* @brief The time stepping starts here.
	*/
	/** If the starting time is not zero, please setup the restart time step or read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		water_block.updateCellLinkedList();
		water_block_complex.updateConfiguration();
	}

	/** Output the start states of bodies. */
	write_water_block_states.writeToFile(0);
	/** Output the mechanical energy of fluid. */
	write_water_mechanical_energy.writeToFile(0);

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 20.0;
	//time step size for output file
	Real D_Time = End_Time / 20.0;
	Real dt = 0.0; //default acoustic time step sizes

	//output for initial particles, global data
	write_water_block_states.writeToFile();

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//computation loop starts
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		//integrate time (loop) until the next output time
		while (integration_time < D_Time)
		{

			//acceleration due to viscous force and gravity
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();

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

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			fluid_observer_contact.updateConfiguration();
			write_recorded_water_pressure.writeToFile(number_of_iterations);
		}

		write_water_mechanical_energy.writeToFile(number_of_iterations);

		tick_count t2 = tick_count::now();
		write_water_block_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	/*
	write_water_mechanical_energy.newResultTest();
	write_recorded_water_pressure.newResultTest();
	*/

	return 0;
}
