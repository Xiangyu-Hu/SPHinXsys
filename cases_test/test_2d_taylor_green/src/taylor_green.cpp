/**
 * @file 	taylor_green.cpp
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
Real DL = 1.0; 						/**< box length. */
Real DH = 1.0; 						/**< box height. */

Real particle_spacing_ref = 1.0/100.0; 		/**< Initial reference particle spacing. */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real U_f = 1.0;							/**< Characteristic velocity. */
Real c_f = 10.0*U_f;					/**< Reference sound speed. */
Real Re = 100;							/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re;		/**< Dynamics viscosity. */

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
 * application dependent initial condition 
 */
class TaylorGreenInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	TaylorGreenInitialCondition(FluidBody *water)
		: fluid_dynamics::FluidInitialCondition(water) {};
protected:
	void Update(size_t index_i, Real dt) override 
	{
		/** initial velocity profile */
		vel_n_[index_i][0] = -cos(2.0 * Pi * pos_n_[index_i][0]) *
				sin(2.0 * Pi * pos_n_[index_i][1]);
		vel_n_[index_i][1] = sin(2.0 * Pi * pos_n_[index_i][0]) *
				cos(2.0 * Pi * pos_n_[index_i][1]);
	}
};
/**
 * @brief 	Main program starts here.
 */
int main(int ac, char* av[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem system(Vec2d(0), Vec2d(DL, DH), particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = false;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);
	
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", 0);
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/** topology */
	SPHBodyComplexRelation* water_block_complex = new SPHBodyComplexRelation(water_block, {});
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /**
	  * @brief 	Methods used only once.
	  */

	/** Initial velocity field */
	TaylorGreenInitialCondition setup_taylor_green_velocity(water_block);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block);
	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding_x(water_block, 0);
	/** Periodic bounding in y direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding_y(water_block, 1);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition_x(water_block, 0);
	/** Periodic BCs in y direction. */
	PeriodicConditionInAxisDirection 	periodic_condition_y(water_block, 1);
	
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_density(water_block_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize 	get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. 
	  * The other reason is that we are using transport velocity formulation, 
	  * which will also introduce numerical disspation slightly. */
	fluid_dynamics::PressureRelaxationFirstHalf pressure_relaxation_first_half(water_block_complex);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann pressure_relaxation_second_half(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAcceleration 	viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityFormulation transport_velocity_formulation(water_block_complex);
	/**
	 * @brief Output.
	 */
	In_Output in_output(system);
	/** Output the body states. */
	WriteBodyStatesToVtu 	write_body_states(in_output, system.real_bodies_);
	/** Write the particle reload files. */
	WriteReloadParticle 		write_particle_reload_files(in_output, { water_block });
	/** Output the body states for restart simulation. */
	ReadRestart				read_restart_files(in_output, system.real_bodies_);
	WriteRestart			write_restart_files(in_output, system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	WriteTotalMechanicalEnergy 	write_total_mechanical_energy(in_output, water_block, new Gravity(Vec2d(0)));
	/** Output the maximum speed of the fluid body. */
	WriteMaximumSpeed write_maximum_speed(in_output, water_block);
	/**
	 * @brief Setup geomtry and initial conditions
	 */
	 /** Using relaxed particle distribution if needed. */
	if (system.reload_particles_) {
		unique_ptr<ReadReloadParticle>
			reload_insert_body_particles(new ReadReloadParticle(in_output, { water_block }, { water_block->GetBodyName() }));
		reload_insert_body_particles->ReadFromFile();
	}
	setup_taylor_green_velocity.exec();
	system.initializeSystemCellLinkedLists();
	periodic_condition_x.parallel_exec();
	periodic_condition_y.parallel_exec();
	system.initializeSystemConfigurations();
	/**
	 * @brief The time stepping starts here.
	 */
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(system.restart_step_);
		water_block->updateCellLinkedList();
		periodic_condition_x.parallel_exec();
		periodic_condition_y.parallel_exec();
		water_block_complex->updateConfiguration();
	}
	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Output the mechanical energy of fluid. */
	write_total_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	 * @brief 	Basic parameters.
	 */
	int number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval*10;
	Real End_Time = 5.0; 	/**< End time. */
	Real D_Time = 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * @brief 	Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force. */
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);
			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				relaxation_time += dt;
				integration_time += dt;
				pressure_relaxation_first_half.parallel_exec(dt);
				pressure_relaxation_second_half.parallel_exec(dt);
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0) {
					write_restart_files.WriteToFile(Real(number_of_iterations));
				}
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_bounding_x.parallel_exec();
			periodic_bounding_y.parallel_exec();
			water_block->updateCellLinkedList();
			periodic_condition_x.parallel_exec();
			periodic_condition_y.parallel_exec();
			water_block_complex->updateConfiguration();
		}

		tick_count t2 = tick_count::now();
		write_total_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
		write_maximum_speed.WriteToFile(GlobalStaticVariables::physical_time_);
		write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << endl;

	write_particle_reload_files.WriteToFile();
	return 0;
}
