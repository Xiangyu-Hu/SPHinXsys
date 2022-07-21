/**
 * @file 	lees_edwards.cpp
 * @brief 	2D lees_edwards boundary condition example.
 * @details test lees_edwards boundary condition.
 * @author 	
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
Real DL = 1.0;					   /**< box length. */
Real DH = 1.0;					   /**< box height. */
Real resolution_ref = 1.0 / 100.0; /**< Global reference resolution. */
Real sr = 0.5;             /**< shear rate. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0), Vec2d(DL, DH));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;					/**< Reference density of fluid. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Reference sound speed. */
Real Re = 100;						/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */


Vec2d ball_center_1(0.5, 0.5);
Real ball_radius = 0.1;
Real gravity_g = 0.0;
Real rho0_s = 1.0;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.0;
//Real physical_viscosity = 10000.0;

// observer location
StdVec<Vecd> observation_location_1 = {ball_center_1};
/**
 * @brief 	Fluid body shape definition.
 */
 /** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	// geometry
	std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
	water_block_shape.push_back(Vecd(0.0, DH));
	water_block_shape.push_back(Vecd(DL, DH));
	water_block_shape.push_back(Vecd(DL, 0.0));
	water_block_shape.push_back(Vecd(0.0, 0.0));

	return water_block_shape;
}
/**
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. 
		std::vector<Vecd> water_block_shape;
		water_block_shape.push_back(Vecd(0.0, 0.0));
		water_block_shape.push_back(Vecd(0.0, DH));
		water_block_shape.push_back(Vecd(DL, DH));
		water_block_shape.push_back(Vecd(DL, 0.0));
		water_block_shape.push_back(Vecd(0.0, 0.0));
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
*/
/** Water block shape definition */
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon.addACircle(ball_center_1, ball_radius, 100, ShapeBooleanOps::sub);
		add<MultiPolygonShape>(multi_polygon);
	}
};

/**
 * apply initial uniform shear velocity profile 
 */
class ShearFlowInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit ShearFlowInitialCondition(FluidBody &water)
		: fluid_dynamics::FluidInitialCondition(water){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		/** initial velocity profile */
		vel_n_[index_i][0] = sr*(pos_n_[index_i][1]-0.5*DH);
	}
};
class FreeBall : public MultiPolygonShape
{
public:
	explicit FreeBall(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addACircle(ball_center_1, ball_radius, 100, ShapeBooleanOps::add);
	}
};

/**
 * @brief 	Main program starts here.
 */
int main(int ac, char *av[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
		/** Tag for running particle relaxation for the initially body-fitted distribution */
	sph_system.run_particle_relaxation_ = false;
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition */
	sph_system.restart_step_ = 0;
	//handle command line arguments
	sph_system.handleCommandlineOptions(ac, av);
	/** output environment. */
	InOutput in_output(sph_system);
	/**
	 * @brief create body, particle and material property.
	 */
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	// Using relaxed particle distribution if needed
	sph_system.reload_particles_
		? water_block.generateParticles<ParticleGeneratorReload>(in_output, water_block.getBodyName())
		: water_block.generateParticles<ParticleGeneratorLattice>();
	/** topology */
	//BodyRelationInner water_block_inner(water_block);
	
	SolidBody free_ball(sph_system, makeShared<FreeBall>("FreeBall"));
	free_ball.defineBodyLevelSetShape();
	free_ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	(!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
		? free_ball.generateParticles<ParticleGeneratorReload>(in_output, free_ball.getBodyName())
		: free_ball.generateParticles<ParticleGeneratorLattice>();

	ObserverBody free_ball_observer(sph_system, "FreeBallObserver");
	free_ball_observer.generateParticles<ObserverParticleGenerator>(observation_location_1);
	
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	/** Initial velocity field */
	ShearFlowInitialCondition initial_condition(water_block);
		//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner water_block_inner(water_block);
	ComplexBodyRelation water_block_complex(water_block_inner, {&free_ball});
	BodyRelationContact free_ball_contact(free_ball, {&water_block});
    BodyRelationInner free_ball_inner(free_ball);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		BodyRelationInner free_ball_inner(free_ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation.
		//----------------------------------------------------------------------
		RandomizeParticlePosition free_ball_random_particles(free_ball);
		relax_dynamics::RelaxationStepInner free_ball_relaxation_step_inner(free_ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		BodyStatesRecordingToVtp write_ball_state(in_output, sph_system.real_bodies_);
		ReloadParticleIO write_particle_reload_files(in_output, {&free_ball});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		free_ball_random_particles.parallel_exec(0.25);
		write_ball_state.writeToFile(0);
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begins.
		//----------------------------------------------------------------------
		int ite = 0;
		int relax_step = 1000;
		while (ite < relax_step)
		{
			free_ball_relaxation_step_inner.exec();
			ite += 1;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_ball_state.writeToFile(ite);
			}
		}
		std::cout << "The physics relaxation process of ball particle finish !" << std::endl;
		write_particle_reload_files.writeToFile(0);
		return 0;
	}

	SimpleDynamics<NormalDirectionFromBodyShape> free_ball_normal_direction(free_ball);
	Gravity gravity(Vecd(0.0, 0.0));
	TimeStepInitialization free_ball_initialize_timestep(free_ball, gravity);
	solid_dynamics::CorrectConfiguration free_ball_corrected_configuration(free_ball_inner);
	solid_dynamics::AcousticTimeStepSize free_ball_get_time_step_size(free_ball);
	/** stress relaxation for the balls. */
	solid_dynamics::StressRelaxationFirstHalf free_ball_stress_relaxation_first_half(free_ball_inner);
	solid_dynamics::StressRelaxationSecondHalf free_ball_stress_relaxation_second_half(free_ball_inner);
	/** Algorithms for solid-solid contact. */
	//solid_dynamics::ContactDensitySummation free_ball_update_contact_density(free_ball_contact);
	//solid_dynamics::ContactForce free_ball_compute_solid_contact_forces(free_ball_contact);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	/** Initialize particle acceleration. */
	TimeStepInitialization time_step_initialization(water_block);
	/** Lees Edwards BCs in y direction */
    LeesEdwardsConditionInAxisDirectionUsingGhostParticles periodic_condition_y(water_block, yAxis);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition_x(water_block, xAxis);

	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationInner update_density_by_summation(water_block_inner);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. 
	  * The other reason is that we are using transport velocity formulation, 
	  * which will also introduce numerical dissipation slightly. */
	fluid_dynamics::PressureRelaxationInner pressure_relaxation(water_block_inner);
	/** Pressure relaxation using Verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation_wall(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation_wall(water_block_complex);

    //pressure_relaxation.pre_processes_.push_back(&periodic_condition_x.ghost_update_);
    pressure_relaxation.pre_processes_.push_back(&periodic_condition_y.ghost_update_);
	fluid_dynamics::DensityRelaxationRiemannInner density_relaxation(water_block_inner);
    //density_relaxation.pre_processes_.push_back(&periodic_condition_x.ghost_update_);
    density_relaxation.pre_processes_.push_back(&periodic_condition_y.ghost_update_);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationInner viscous_acceleration(water_block_inner);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionInner transport_velocity_correction(water_block_inner);

	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	/** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_inserted_body(free_ball_contact);
	/** Computing viscous force acting on wall with wall model. */
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_inserted_body(free_ball_contact);

	/**
	 * @brief Output.
	 */
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	/** Write the particle reload files. */
	ReloadParticleIO write_particle_reload_files(in_output, {&water_block});
	/** Output the body states for restart simulation. */
	RestartIO restart_io(in_output, sph_system.real_bodies_);
	/** Output the maximum speed of the fluid body. */
	RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<MaximumSpeed>>
		write_maximum_speed(in_output, water_block);
	RegressionTestTimeAveraged<BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid>>
		write_total_viscous_force_on_inserted_body(in_output, free_ball);
	BodyReducedQuantityRecording<solid_dynamics::TotalViscousForceOnSolid>
		write_total_force_on_inserted_body(in_output, free_ball);
	/**
	 * @brief Setup geometry and initial conditions
	 */
	initial_condition.exec();
	sph_system.initializeSystemCellLinkedLists();
    // initial periodic boundary condition
    //periodic_condition_x.ghost_creation_.parallel_exec();
    periodic_condition_y.ghost_creation_.parallel_exec();
	periodic_condition_x.update_cell_linked_list_.parallel_exec();
	//periodic_condition_y.update_cell_linked_list_.parallel_exec();
	sph_system.initializeSystemConfigurations();
	free_ball_corrected_configuration.parallel_exec();
	/** initialize surface normal direction for the insert body. */
	free_ball_normal_direction.parallel_exec();
	/**
	 * @brief The time stepping starts here.
	 */
	/** If the starting time is not zero, please setup the restart time step or read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		free_ball.updateCellLinkedList();
		water_block.updateCellLinkedList();
        //periodic_condition_x.ghost_creation_.parallel_exec();
        periodic_condition_y.ghost_creation_.parallel_exec();
		periodic_condition_x.update_cell_linked_list_.parallel_exec();
		//periodic_condition_y.update_cell_linked_list_.parallel_exec();
		water_block_inner.updateConfiguration();
		free_ball_contact.updateConfiguration();
	}
	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 20.0; /**< End time. */
	Real D_Time = 0.1;	 /**< Time stamps for output of body states. */
	Real dt = 0.0;		 /**< Default acoustic time step sizes. */
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
			time_step_initialization.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);
			/** FSI for viscous force. */
			fluid_viscous_force_on_inserted_body.parallel_exec();

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				free_ball_initialize_timestep.parallel_exec();
				//free_ball_update_contact_density.parallel_exec();
				free_ball_stress_relaxation_first_half.parallel_exec(dt);
				free_ball_stress_relaxation_second_half.parallel_exec(dt);
				free_ball.updateCellLinkedList();
				Real dt_free = free_ball_get_time_step_size.parallel_exec();

				//avoid possible smaller acoustic time step size for viscous flow
				Real dt_fluid = SMIN(get_fluid_time_step_size.parallel_exec(), Dt);
                dt = SMIN(dt_fluid,dt_free);

				relaxation_time += dt;
				integration_time += dt;

				pressure_relaxation.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_inserted_body.parallel_exec();
				density_relaxation.parallel_exec(dt);
				GlobalStaticVariables::physical_time_ += dt;

			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
				{
					restart_io.writeToFile(number_of_iterations);
				}
			}
			number_of_iterations++;
             
			/** Water block configuration and periodic condition. */
		    periodic_condition_y.bounding_.parallel_exec();
			periodic_condition_x.bounding_.parallel_exec();
			water_block.updateCellLinkedList();
            //periodic_condition_x.ghost_creation_.parallel_exec();
            periodic_condition_y.ghost_creation_.parallel_exec();
			periodic_condition_x.update_cell_linked_list_.parallel_exec();
			//periodic_condition_y.update_cell_linked_list_.parallel_exec();
			water_block_inner.updateConfiguration();
			free_ball_contact.updateConfiguration();
		}

		tick_count t2 = tick_count::now();
		write_maximum_speed.writeToFile(number_of_iterations);
		body_states_recording.writeToFile();
		write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
		write_total_force_on_inserted_body.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;

	write_particle_reload_files.writeToFile();
	if (!sph_system.reload_particles_)
	{
		write_maximum_speed.newResultTest();
	}

	return 0;
}
