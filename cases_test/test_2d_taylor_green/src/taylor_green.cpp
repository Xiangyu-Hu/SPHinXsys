/**
 * @file 	Taylor_Green vortex flow.cpp
 * @brief 	2D taylor_green vortex flow example.
 * @details This is the one of the basic test cases.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			this code developed based one the beta version
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
Real DL = 1.0; 						/**< Tank length. */
Real DH = 1.0; 						/**< Tank height. */

Real particle_spacing_ref = 1.0/100; 		/**< Initial reference particle spacing. */
Real BW = 4.0 * particle_spacing_ref;
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real gravity_g = 0.0;					/**< Gravity force of fluid. */
Real U_f = 1.0;							/**< Characteristic velocity. */
Real c_f = 10.0*U_f;					/**< Reference sound speed. */
Real Re = 100;
Real mu_f = rho0_f * U_f * DL / Re;		/**< Dynamics visocisty. */
Real k_f = 0.0;

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
 * application dependent initial condition 
 */
class TaylorGreenInitialCondition
	: public fluid_dynamics::WeaklyCompressibleFluidInitialCondition
{
public:
	TaylorGreenInitialCondition(FluidBody *water)
		: fluid_dynamics::WeaklyCompressibleFluidInitialCondition(water) {};
protected:
	void Update(size_t index_particle_i, Real dt) override 
	{
		/** first set all particle at rest*/
		fluid_dynamics::WeaklyCompressibleFluidInitialCondition::Update(index_particle_i, dt);

		/** initial velocity profile */
		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		FluidParticleData &fluid_data_i = particles_->fluid_particle_data_[index_particle_i];

		base_particle_data_i.vel_n_[0] = -cos(2.0 * pi * base_particle_data_i.pos_n_[0]) * 
				sin(2.0 * pi * base_particle_data_i.pos_n_[1]);
		base_particle_data_i.vel_n_[1] = sin(2.0 * pi * base_particle_data_i.pos_n_[0]) * 
				cos(2.0 * pi * base_particle_data_i.pos_n_[1]);
	};
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
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", 
		0, ParticlesGeneratorOps::lattice);
	WeaklyCompressibleFluid 		fluid("Water", water_block, rho0_f, c_f, mu_f, k_f);
	FluidParticles 	fluid_particles(water_block);
	/**
	 * @brief 	Body contact map.
	 * @details The contact map gives the data conntections between the bodies.
	 * 			Basically the the rang of bidies to build neighbor particle lists.
	 */
	SPHBodyTopology 	body_topology = { { water_block, { } } };
	system.SetBodyTopology(&body_topology);
	/**
	 * @brief 	Simulation set up.
	 */
	system.SetupSPHSimulation();
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Define external force. */
	Gravity 							gravity(Vecd(0.0, 0.0));
	 /**
	  * @brief 	Methods used only once.
	  */

	/** Obtain the initial number density of fluid. */
	fluid_dynamics::InitialNumberDensity 		fluid_initial_number_density(water_block, { });
	/** Reset velocity field */
	TaylorGreenInitialCondition setup_taylor_green_velocity(water_block);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeOtherAccelerations 	initialize_fluid_acceleration(water_block, &gravity);
	/** Periodic bounding. */
	PeriodicBoundingInXDirection 	periodic_bounding_x(water_block);
	PeriodicBoundingInYDirection 	periodic_bounding_y(water_block);
	/** Periodic BCs. */
	PeriodicConditionInXDirection 	periodic_condition_x(water_block);
	PeriodicConditionInYDirection 	periodic_condition_y(water_block);
	
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_desnity(water_block, { });
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize 	get_fluid_adevction_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize 	get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	fluid_dynamics::PressureRelaxationVerlet 	pressure_relaxation(water_block, { }, &gravity);
	/** Computing viscous acceleration. */
	fluid_dynamics::ComputingViscousAcceleration 	viscous_acceleration(water_block, {  });
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrection 	transport_velocity_correction(water_block, { });
	/** Computing vorticity in the flow. */
	/**
	 * @brief 	Methods used for updating data structure.
	 */
	 /** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList			update_cell_linked_list(water_block);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsConfiguration 			update_particle_configuration(water_block);
	/**
	 * @brief Output.
	 */
	In_Output in_output(system);
	/** Output the body states. */
	WriteBodyStatesToVtu 	write_body_states(in_output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	ReadRestart				read_restart_files(in_output, system.real_bodies_);
	WriteRestart			write_restart_files(in_output, system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	WriteWaterMechanicalEnergy 	write_water_mechanical_energy(in_output, water_block, &gravity);

	/**
	 * @brief Setup goematrics and initial conditions
	 */
	setup_taylor_green_velocity.exec();
	periodic_condition_x.parallel_exec();
	periodic_condition_y.parallel_exec();
	update_particle_configuration.parallel_exec();
	fluid_initial_number_density.exec();
	/**
	 * @brief The time stepping starts here.
	 */
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(system.restart_step_);
		update_cell_linked_list.parallel_exec();
		update_particle_configuration.parallel_exec();
	}
	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	 * @brief 	Basic parameters.
	 */
	int number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval*10;
	Real End_Time = 4.0; 	/**< End time. */
	Real D_Time = 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default accoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * @brief 	Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integeral_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			initialize_fluid_acceleration.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					write_restart_files.WriteToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** Water block confifuration and periodic constion. */
			periodic_bounding_x.parallel_exec();
			periodic_bounding_y.parallel_exec();
			update_cell_linked_list.parallel_exec();
			periodic_condition_x.parallel_exec();
			periodic_condition_y.parallel_exec();
			update_particle_configuration.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
		write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << endl;

	return 0;
}
