/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/
/**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
Real resolution_ref = 0.05; //particle spacing
Real BW = resolution_ref * 4; //boundary width
Real DL = 5.366; 						//tank length
Real DH = 2.0; 							//tank height
Real DW = 0.5;							//tank width
Real LL = 2.0; 							//liquid length
Real LH = 1.0; 							//liquid height
Real LW = 0.5; 						//liquid width
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));

//for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0*sqrt(gravity_g * LH);
Real c_f = 10.0*U_f;

/* resolution which control the quality of polygonalmesh created by geometry system */
int resolution(50);

//define the fluid body
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, std::string body_name)
		: FluidBody(system, body_name)
	{
		Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
		Vecd translation_water = halfsize_water;
		mesh_.reset(new ComplexShapeTriangleMesh());
		mesh_->addBrick(halfsize_water, resolution, translation_water, ShapeBooleanOps::add);
		body_shape_ = new ComplexShape(mesh_.get());
	}

private:
	std::unique_ptr<ComplexShapeTriangleMesh> mesh_;

};
/**
 * @brief 	Case dependent material properties definition.
 */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		rho0_ = rho0_f;
		c0_ = c_f;

		assignDerivedMaterialParameters();
	}
};
//define the static solid wall boundary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
		Vecd translation_wall(0.5 * DL, 0.5 * DH, 0.5 * DW);
		Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
		
		mesh_.reset(new ComplexShapeTriangleMesh());
		mesh_->addBrick(halfsize_outer, resolution, translation_wall, ShapeBooleanOps::add);
		mesh_->addBrick(halfsize_inner, resolution, translation_wall, ShapeBooleanOps::sub);
		body_shape_ = new ComplexShape(mesh_.get());
	}
private:
	std::unique_ptr<ComplexShapeTriangleMesh> mesh_;
};

//define an observer body
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem &system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		//add observation point
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.01, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.1, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.2, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.24, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.252, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.266, 0.5 * DW), 0.0));
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
	WaterBlock *water_block = new WaterBlock(system, "WaterBody");
	//Configuration of Materials
	WaterMaterial *water_material = new WaterMaterial();
	//creat fluid particles
	FluidParticles fluid_particles(water_block, water_material);

	//the wall boundary
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall");
	//creat solid particles 
	SolidParticles wall_particles(wall_boundary);

	FluidObserver *fluid_observer = new FluidObserver(system, "Fluidobserver");
	//create observer particles 
	BaseParticles observer_particles(fluid_observer);

	//define external force
	Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));


	/** topology */
	ComplexBodyRelation* water_block_complex = new ComplexBodyRelation(water_block, { wall_boundary });
	BodyRelationContact* fluid_observer_contact = new BodyRelationContact(fluid_observer, { water_block });

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------
	//-------- common particle dynamics ----------------------------------------
	TimeStepInitialization 	initialize_a_fluid_step(water_block, &gravity);
	//-------- fluid dynamics --------------------------------------------------
	//evaluation of density by summation approach
	fluid_dynamics::DensitySummationFreeSurfaceComplex	update_density_by_summation(water_block_complex);
	//time step size without considering sound wave speed
	fluid_dynamics::AdvectionTimeStepSize	get_fluid_advection_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(water_block);

	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationRiemannWithWall	pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall	density_relaxation(water_block_complex);

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_water_block_states(in_output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO		restart_io(in_output, system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	BodyReducedQuantityRecording<TotalMechanicalEnergy> 	
		write_water_mechanical_energy(in_output, water_block, &gravity);
	/** output the observed data from fluid body. */
	ObservedQuantityRecording<indexScalar, Real> write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);
	//-------------------------------------------------------------------
	//from here the time stepping begines
	//-------------------------------------------------------------------

	/**
	 * @brief Setup geemetrics and initial conditions
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	/**
	* @brief The time stepping starts here.
	*/
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		water_block->updateCellLinkedList();
		water_block_complex->updateConfiguration();
	}
	
	/** Output the start states of bodies. */
	write_water_block_states.writeToFile(0);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.writeToFile(0);

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 20.0;
	//time step size for ouput file
	Real D_Time = End_Time/20.0;
	Real Dt = 0.0;//default advection time step sizes
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
			Dt = get_fluid_advection_time_step_size.parallel_exec();
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

			water_block->updateCellLinkedList();
			water_block_complex->updateConfiguration();
			fluid_observer_contact->updateConfiguration();
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

	return 0;
}
