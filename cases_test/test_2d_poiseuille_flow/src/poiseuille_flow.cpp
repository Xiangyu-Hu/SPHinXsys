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
Real mu_f = 1.0e-6;							/**< Visocitiy. */
Real U_f = gravity_g * DH * DH / mu_f;		/**< Characteristic velocity. */
Real c_f = 10.0*U_f;						/**< Reference sound speed. */
/**
 * @brief 	Fluid body definition.
 */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: FluidBody(system, body_name, refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> water_block_shape;
		water_block_shape.push_back(Point(0.0, 0.0));
		water_block_shape.push_back(Point(0.0, DH));
		water_block_shape.push_back(Point(DL, DH));
		water_block_shape.push_back(Point(DL, 0.0));
		water_block_shape.push_back(Point(0.0, 0.0));
		Geometry *water_block_geometry = new Geometry(water_block_shape);
		body_region_.add_geometry(water_block_geometry, RegionBooleanOps::add);

		body_region_.done_modeling();
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
	WallBoundary(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> outer_wall_shape;
		outer_wall_shape.push_back(Point(-BW, -BW));
		outer_wall_shape.push_back(Point(-BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, -BW));
		outer_wall_shape.push_back(Point(-BW, -BW));
		Geometry *outer_wall_geometry = new Geometry(outer_wall_shape);
		body_region_.add_geometry(outer_wall_geometry, RegionBooleanOps::add);

		std::vector<Point> inner_wall_shape;
		inner_wall_shape.push_back(Point(-2.0 * BW, 0.0));
		inner_wall_shape.push_back(Point(-2.0 * BW, DH));
		inner_wall_shape.push_back(Point(DL + 2.0 * BW, DH));
		inner_wall_shape.push_back(Point(DL + 2.0 * BW, 0.0));
		inner_wall_shape.push_back(Point(-2.0 * BW, 0.0));
		Geometry *inner_wall_geometry = new Geometry(inner_wall_shape);
		body_region_.add_geometry(inner_wall_geometry, RegionBooleanOps::sub);

		body_region_.done_modeling();
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
	WaterBlock *water_block 
		= new WaterBlock(system, "WaterBody", 0, ParticlesGeneratorOps::lattice);
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary 
		= new WallBoundary(system, "Wall",	0, ParticlesGeneratorOps::lattice);
	SolidParticles 					solid_particles(wall_boundary);
	/**
	 * @brief 	Body contact map.
	 * @details The contact map gives the data conntections between the bodies.
	 * 			Basically the the rang of bidies to build neighbor particle lists.
	 */
	SPHBodyTopology 	body_topology = { { water_block, { wall_boundary } },
										  { wall_boundary, {} }};
	system.SetBodyTopology(&body_topology);

	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Define external force. */
	Gravity gravity(Vecd(gravity_g, 0.0));
	 /**
	  * @brief 	Methods used only once.
	  */
	/** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_boundary, {});
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block, &gravity);
	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(water_block, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(water_block, 0);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 		update_fluid_density(water_block, { wall_boundary });
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize 	get_fluid_adevction_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm without Riemann solver for viscous flows. */
	fluid_dynamics::PressureRelaxationFirstHalf 
		pressure_relaxation_first_half(water_block, { wall_boundary });
	/** Pressure relaxation algorithm by using position verlet time stepping. */
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(water_block, { wall_boundary });
	/** Computing viscous acceleration. */
	fluid_dynamics::ComputingViscousAcceleration 	
		viscous_acceleration(water_block, { wall_boundary});
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrection transport_velocity_correction(water_block, { wall_boundary});
	/** Computing implicit viscous acceleration. */
	fluid_dynamics::ImplicitComputingViscousAcceleration 	
		implicit_viscous_acceleration(water_block, { wall_boundary});
	/**
	 * @brief 	Methods used for updating data structure.
	 */
	 /** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList	update_cell_linked_list(water_block);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsConfiguration 	update_particle_configuration(water_block);
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
	 * @brief Setup goematrics and initial conditions.
	 */
	system.InitializeSystemCellLinkedLists();
	system.InitializeSystemConfigurations();
	get_wall_normal.exec();
	/**
	 * @brief The time stepping starts here.
	 */
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(system.restart_step_);
		periodic_condition.parallel_exec();
		update_cell_linked_list.parallel_exec();
		update_particle_configuration.parallel_exec();
	}
	/** Pre-simulation*/
	periodic_condition.parallel_exec();
	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	 * @brief 	Basic parameters.
	 */
	int number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval*10;
	Real End_Time 		= 20.0; 	/**< End time. */
	Real Output_Time 	= 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default accoustic time step sizes. */
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
		Real integeral_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integeral_time < Output_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			//viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);
			interval_computing_time_step += tick_count::now() - time_instance;
			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation_first_half.parallel_exec(dt);
				implicit_viscous_acceleration.parallel_exec(dt);
				pressure_relaxation_second_half.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
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
			periodic_bounding.parallel_exec();
			update_cell_linked_list.parallel_exec();
			periodic_condition.parallel_exec();
			update_particle_configuration.parallel_exec();
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
