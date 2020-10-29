/**
 * @file 	poiseuille_flow.cpp
 * @brief 	2D poiseuille flow example
 * @details This is the one of the basic test cases for validating the splitting-implicit SPH method.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
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
Real DL = 1.0e-3; 						/**< Tank length. */
Real DH = 1.0e-3; 						/**< Tank height. */
Real particle_spacing_ref = DH / 20.0; 		/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1000.0;						/**< Reference density of fluid. */
Real gravity_g = 1.0e-4;					/**< Gravity force of fluid. */
Real mu_f = 1.0e-6;							/**< Viscosity. */
Real U_f = gravity_g * DH * DH / mu_f;		/**< Characteristic velocity. */
Real c_f = 10.0*U_f;						/**< Reference sound speed. */
/**
 * @brief 	Fluid body definition.
 */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, string body_name,	int refinement_level)
		: FluidBody(system, body_name, refinement_level)
	{
		/** Geomtry definition. */
		std::vector<Point> water_block_shape;
		water_block_shape.push_back(Point(0.0, 0.0));
		water_block_shape.push_back(Point(0.0, DH));
		water_block_shape.push_back(Point(DL, DH));
		water_block_shape.push_back(Point(DL, 0.0));
		water_block_shape.push_back(Point(0.0, 0.0));
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
	WaterMaterial()	: WeaklyCompressibleFluid()
	{
		rho_0_ = rho0_f;
		c_0_ = c_f;
		mu_ = mu_f;

		assignDerivedMaterialParameters();
	}
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geomtry definition. */
		std::vector<Point> outer_wall_shape;
		outer_wall_shape.push_back(Point(-BW, -BW));
		outer_wall_shape.push_back(Point(-BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, -BW));
		outer_wall_shape.push_back(Point(-BW, -BW));
		std::vector<Point> inner_wall_shape;
		inner_wall_shape.push_back(Point(-2.0 * BW, 0.0));
		inner_wall_shape.push_back(Point(-2.0 * BW, DH));
		inner_wall_shape.push_back(Point(DL + 2.0 * BW, DH));
		inner_wall_shape.push_back(Point(DL + 2.0 * BW, 0.0));
		inner_wall_shape.push_back(Point(-2.0 * BW, 0.0));
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
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
	SPHSystem system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", 0);
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall",	0);
	SolidParticles 	wall_particles(wall_boundary);
	/** topology */
	SPHBodyComplexRelation* water_block_complex = new SPHBodyComplexRelation(water_block, { wall_boundary });
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Define external force. */
	Gravity gravity(Vecd(gravity_g, 0.0));
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block, &gravity);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList 	periodic_condition(water_block, 0);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 		update_fluid_density(water_block_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm without Riemann solver for viscous flows. */
	fluid_dynamics::PressureRelaxationFirstHalf 
		pressure_relaxation_first_half(water_block_complex);
	/** Pressure relaxation algorithm by using position verlet time stepping. */
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAcceleration 	
		viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityFormulation transport_velocity_formulation(water_block_complex);
	/**
	 * @brief Output.
	 */
	In_Output in_output(system);
	/** Output the body states. */
	WriteBodyStatesToVtu write_body_states(in_output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	ReadRestart		read_restart_files(in_output, system.real_bodies_);
	WriteRestart	write_restart_files(in_output, system.real_bodies_);
	/**
	 * @brief Setup geomtry and initial conditions.
	 */
	system.initializeSystemCellLinkedLists();
	periodic_condition.update_cell_linked_list_.parallel_exec();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	/**
	 * @brief The time stepping starts here.
	 */
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(system.restart_step_);
		water_block->updateCellLinkedList();
		periodic_condition.update_cell_linked_list_.parallel_exec();
		water_block_complex->updateConfiguration();
	}
	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval*10;
	Real End_Time 		= 20.0; 	/**< End time. */
	Real Output_Time 	= 0.1;		/**< Time stamps for output of body states. */
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
		while (integration_time < Output_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			//viscous_acceleration.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);
			interval_computing_time_step += tick_count::now() - time_instance;
			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation_first_half.parallel_exec(dt);
				viscous_acceleration.parallel_exec(dt);
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
			/** Water block configuration and periodic condition. */
			periodic_condition.bounding_.parallel_exec();
			water_block->updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
			water_block_complex->updateConfiguration();
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
