/**
 * @file 	Dambreak.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "sphinxsys.h"	//SPHinXsys Library.

using namespace SPH;	//Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366; 						/**< Tank length. */
Real DH = 5.366; 						/**< Tank height. */
Real LL = 2.0; 							/**< Liquid colume length. */
Real LH = 1.0; 							/**< Liquid colume height. */
Real particle_spacing_ref = 0.025; 		/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real gravity_g = 1.0;					/**< Gravity force of fluid. */
Real U_max = 2.0*sqrt(gravity_g*LH);	/**< Characteristic velocity. */
Real c_f = 10.0* U_max;					/**< Reference sound speed. */
/** create a water block shape */
std::vector<Vecd> CreatWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.0, 0.0));
	water_block_shape.push_back(Vecd(0.0, LH));
	water_block_shape.push_back(Vecd(LL, LH));
	water_block_shape.push_back(Vecd(LL, 0.0));
	water_block_shape.push_back(Vecd(0.0, 0.0));
	return water_block_shape;
}
/** create outer wall shape */
std::vector<Vecd> CreatOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));

	return outer_wall_shape;
}
/** create inner wall shape */
std::vector<Vecd> CreatInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

	return inner_wall_shape;
}
//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& sph_system, std::string body_name)
		: FluidBody(sph_system, body_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> water_block_shape = CreatWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Case dependent material properties definition.
//----------------------------------------------------------------------
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		/** Basic material parameters*/
		rho0_ = rho0_f;
		c0_ = c_f;

		/** Compute the derived material parameters*/
		assignDerivedMaterialParameters();
	}
};
//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, std::string body_name)
		: SolidBody(sph_system, body_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> outer_shape = CreatOuterWallShape();
		std::vector<Vecd> inner_shape = CreatInnerWallShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_shape, ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Fluid observer body definition.
//----------------------------------------------------------------------
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem &sph_system, std::string body_name)
		: FictitiousBody(sph_system, body_name)
	{
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.2), 0.0));
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	GlobalStaticVariables::physical_time_ = 0.0;
	sph_system.handleCommandlineOptions(ac, av);
	In_Output in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	WaterBlock		*water_block = new WaterBlock(sph_system, "WaterBody");
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);

	WallBoundary	*wall_boundary = new WallBoundary(sph_system, "Wall");
	SolidParticles	wall_particles(wall_boundary);

	FluidObserver	*fluid_observer = new FluidObserver(sph_system, "Fluidobserver");
	BaseParticles	observer_particles(fluid_observer);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexBodyRelation* water_block_complex = new ComplexBodyRelation(water_block, { wall_boundary });
	ContactBodyRelation* fluid_observer_contact = new ContactBodyRelation(fluid_observer, { water_block });
	//----------------------------------------------------------------------
	//	Define all numerical methods which are used in this case.
	//----------------------------------------------------------------------
	Gravity	gravity(Vecd(0.0, -gravity_g));
	InitializeATimeStep		initialize_a_fluid_step(water_block, &gravity);
	fluid_dynamics::DensitySummationFreeSurfaceComplex 	update_density_by_summation(water_block_complex);
	fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step_size(water_block, U_max);
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	//----------------------------------------------------------------------
	//	out_put body states and observations
	//----------------------------------------------------------------------
	WriteBodyStatesToVtu	write_body_states(in_output, sph_system.real_bodies_);
	RestartIO	restart_io(in_output, sph_system.real_bodies_);
	WriteBodyReducedQuantity<TotalMechanicalEnergy> 	
		write_water_mechanical_energy(in_output, water_block, &gravity);
	WriteAnObservedQuantity<indexScalar, Real>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	initialize the simulation with cell linked list and configuration
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	//----------------------------------------------------------------------
	//	load restart file if necessary
	//----------------------------------------------------------------------
 	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block->updateCellLinkedList();
		water_block_complex->updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	first out put before the simulation
	//----------------------------------------------------------------------
 	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 20.0; 	/**< End time. */
	Real D_Time = 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	statistics for computing CPU time
	//----------------------------------------------------------------------
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
			/** outer loop for dual-time criteria time-stepping. */
			time_instance = tick_count::now();
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			interval_computing_time_step += tick_count::now() - time_instance;

			/** inner loop for dual-time criteria time-stepping.  */
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

			/** screen output, write body reduced values and restart files  */
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % observation_sample_interval == 0) {
					write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
					write_recorded_water_pressure.WriteToFile(GlobalStaticVariables::physical_time_);
				}
				if (number_of_iterations % restart_output_interval == 0)
					restart_io.WriteToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();
			water_block->updateCellLinkedList();
			water_block_complex->updateConfiguration();
			fluid_observer_contact->updateConfiguration();
			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
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
