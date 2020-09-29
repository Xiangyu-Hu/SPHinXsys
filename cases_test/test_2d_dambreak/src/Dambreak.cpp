/**
 * @file 	Dambreak.cpp
 * @brief 	2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 * @version 0.1
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
Real DL = 5.366; 						/**< Tank length. */
Real DH = 5.366; 						/**< Tank height. */
Real LL = 2.0; 							/**< Liquid colume length. */
Real LH = 1.0; 							/**< Liquid colume height. */
Real particle_spacing_ref = 0.025; 		/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real gravity_g = 1.0;					/**< Gravity force of fluid. */
Real U_max = 2.0*sqrt(gravity_g*LH);		/**< Characteristic velocity. */
Real c_f = 10.0* U_max;					/**< Reference sound speed. */
/** create a water block shape */
std::vector<Point> CreatWaterBlockShape()
{
	//geometry
	std::vector<Point> water_block_shape;
	water_block_shape.push_back(Point(0.0, 0.0));
	water_block_shape.push_back(Point(0.0, LH));
	water_block_shape.push_back(Point(LL, LH));
	water_block_shape.push_back(Point(LL, 0.0));
	water_block_shape.push_back(Point(0.0, 0.0));
	return water_block_shape;
}
/** create outer wall shape */
std::vector<Point> CreatOuterWallShape()
{
	std::vector<Point> outer_wall_shape;
	outer_wall_shape.push_back(Point(-BW, -BW));
	outer_wall_shape.push_back(Point(-BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, -BW));
	outer_wall_shape.push_back(Point(-BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Point> CreatInnerWallShape()
{
	std::vector<Point> inner_wall_shape;
	inner_wall_shape.push_back(Point(0.0, 0.0));
	inner_wall_shape.push_back(Point(0.0, DH));
	inner_wall_shape.push_back(Point(DL, DH));
	inner_wall_shape.push_back(Point(DL, 0.0));
	inner_wall_shape.push_back(Point(0.0, 0.0));

	return inner_wall_shape;
}
/**
*@brief 	Fluid body definition.
*/
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& sph_system, string body_name, int refinement_level)
		: FluidBody(sph_system, body_name, refinement_level)
	{
		/** Geomtry definition. */
		std::vector<Point> water_block_shape = CreatWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
/**
 * @brief 	Case dependent material properties definition.
 */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		/** Basic material parameters*/
		rho_0_ = rho0_f;
		c_0_ = c_f;

		/** Compute the derived material parameters*/
		assignDerivedMaterialParameters();
	}
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, string body_name, int refinement_level)
		: SolidBody(sph_system, body_name, refinement_level)
	{
		/** Geomtry definition. */
		std::vector<Point> outer_shape = CreatOuterWallShape();
		std::vector<Point> inner_shape = CreatInnerWallShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_shape, ShapeBooleanOps::sub);
	}
};
/**
 * @brief 	Fluid observer body definition.
 */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem &sph_system, string body_name, int refinement_level)
		: FictitiousBody(sph_system, body_name, refinement_level, 1.3)
	{
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.2), 0.0));
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
	SPHSystem sph_system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	WaterBlock *water_block = new WaterBlock(sph_system, "WaterBody", 0);
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary = new WallBoundary(sph_system, "Wall",	0);
	SolidParticles 					solid_particles(wall_boundary);
	/**
	 * @brief 	Particle and body creation of fluid observer.
	 */
	FluidObserver *fluid_observer = new FluidObserver(sph_system, "Fluidobserver", 0);
	BaseParticles 	observer_particles(fluid_observer);

	/** topology */
	SPHBodyComplexRelation* water_block_complex_relation = new SPHBodyComplexRelation(water_block, { wall_boundary });
	SPHBodyComplexRelation* wall_complex_relation = new SPHBodyComplexRelation(wall_boundary, {});
	SPHBodyContactRelation* fluid_observer_contact_relation = new SPHBodyContactRelation(fluid_observer, { water_block });

	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Define external force. */
	Gravity 							gravity(Vecd(0.0, -gravity_g));
	 /**
	  * @brief 	Methods used only once.
	  */
	/** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_complex_relation);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block, &gravity);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummationFreeSurface 		update_fluid_density(water_block_complex_relation);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 			get_fluid_advection_time_step_size(water_block, U_max);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using position verlet time stepping. */
	fluid_dynamics::PressureRelaxationFirstHalfRiemann 
		pressure_relaxation_first_half(water_block_complex_relation);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann 
		pressure_relaxation_second_half(water_block_complex_relation);

	/**
	 * @brief Output.
	 */
	In_Output in_output(sph_system);
	/** Output the body states. */
	WriteBodyStatesToVtu 		write_body_states(in_output, sph_system.real_bodies_);
	/** Output the body states for restart simulation. */
	ReadRestart		read_restart_files(in_output, sph_system.real_bodies_);
	WriteRestart	write_restart_files(in_output, sph_system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	WriteTotalMechanicalEnergy 	write_water_mechanical_energy(in_output, water_block, &gravity);
	/** output the observed data from fluid body. */
	WriteAnObservedQuantity<Real, FluidParticles, &FluidParticles::p_>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact_relation);

	/** Pre-simulation*/
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	get_wall_normal.exec();

	/**
	 * @brief The time stepping starts here.
	 */
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(sph_system.restart_step_);
		water_block->updateCellLinkedList();
		water_block_complex_relation->updateConfiguration();
	}

	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	 * @brief 	Basic parameters.
	 */
	int number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval*10;
	Real End_Time = 20.0; 	/**< End time. */
	Real D_Time = 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
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
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			interval_computing_time_step += tick_count::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation_first_half.parallel_exec(dt);
				pressure_relaxation_second_half.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

			}
			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					write_restart_files.WriteToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();
			water_block->updateCellLinkedList();
			water_block_complex_relation->updateConfiguration();
			fluid_observer_contact_relation->updateConfiguration();
			interval_updating_configuration += tick_count::now() - time_instance;
		}


		tick_count t2 = tick_count::now();
		write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
		write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		write_recorded_water_pressure.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;

	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << endl;
	cout << fixed << setprecision(9) << "interval_computing_time_step ="
		<< interval_computing_time_step.seconds() << "\n";
	cout << fixed << setprecision(9) << "interval_computing_pressure_relaxation = "
		<< interval_computing_pressure_relaxation.seconds() << "\n";
	cout << fixed << setprecision(9) << "interval_updating_configuration = "
		<< interval_updating_configuration.seconds() << "\n";

	return 0;
}
