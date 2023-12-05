/**
 * @file 	Oscillation drop.cpp
 * @brief 	2D oscillation drop.
 * @author 	Yaru Ren, Chi Zhang and Xiangyu Hu
 */
 /**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
#include "energy.h"
  /**
  * @brief Namespace cite here.
  */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real particle_spacing_ref = 1.0 / 100; 		/**< Initial reference particle spacing.*/
Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1000.0;						/**< Reference density of fluid.*/
Real U_f = 1.0;		                       /**< Characteristic velocity.*/    //2.0*sqrt(gravity_g*LH)
Real c_f = 20.0 * U_f;					  /**< Reference sound speed.*/ //5.0*U_f

//for the circle parameter
Vec2d Circle_center(0, 0);    //circle center
Real  Circle_radius = 1.0; //circle radius 0.1
int   Circle_resolution(100); //resolution

Real DL = 4 * Circle_radius;
Real DH = 4 * Circle_radius;

class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addACircle(Circle_center, Circle_radius, 100, ShapeBooleanOps::add);
	}
};
/**
 * application dependent initial velocity
 */
class InitialVelocity
	: public fluid_dynamics::FluidInitialCondition
{
public:
	InitialVelocity(SPHBody& sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body) {};

	void update(size_t index_particle_i, Real dt)
	{
		/** initial velocity profile */

		vel_[index_particle_i][0] = 1.0 * pos_[index_particle_i][0];
		vel_[index_particle_i][1] = -1.0 * pos_[index_particle_i][1];
		//p_[index_particle_i] = (1.0 - (pos_[index_particle_i][0] * pos_[index_particle_i][0] + pos_[index_particle_i][1] * pos_[index_particle_i][1]));
	}
};
/**
 * application dependent
 */
class ExternalField
	: public fluid_dynamics::FluidInitialCondition
{
public:
	ExternalField(SPHBody& sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body),
		acc_prior_(particles_->acc_prior_) {};
protected:
	void update(size_t index_particle_i, Real dt)
	{
		//acc_prior_[index_particle_i][0] = -1.2 * 1.2 * pos_[index_particle_i][0];
		//acc_prior_[index_particle_i][1] = -1.2 * 1.2 * pos_[index_particle_i][1];
		acc_prior_[index_particle_i][0] = -1.0 * 1.0 * pos_[index_particle_i][0];
		acc_prior_[index_particle_i][1] = -1.0 * 1.0 * pos_[index_particle_i][1];
	}

protected:
	StdLargeVec<Vecd>& acc_prior_;
};

/**
 * @brief 	Main program starts here.
 */
int main(int ac, char* av[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	BoundingBox system_domain_bounds(Vec2d(-DL - BW, -DH - BW), Vec2d(DL + BW, DH + BW));
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	sph_system.handleCommandlineOptions(ac, av);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	sph_system.setRunParticleRelaxation(false);
	/** Tag for computation start with relaxed body fitted particles distribution. */
	sph_system.setReloadParticles(true);
	IOEnvironment io_environment(sph_system);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineBodyLevelSetShape();
	water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	// Using relaxed particle distribution if needed
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
		: water_block.generateParticles<ParticleGeneratorLattice>();

	ObserverBody fluid_observer(sph_system, "FluidObserver");
	StdVec<Vecd> observation_location = { Vecd(0.0, 0.0) };
	fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_block);
	ContactRelation fluid_observer_contact(fluid_observer, { &water_block });

	/** check whether run particle relaxation for body fitted particle distribution. */
	if (sph_system.RunParticleRelaxation())
	{
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		 /** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
		/** Write the particle reload files. */
		ReloadParticleIO write_real_body_particle_reload_files(io_environment, sph_system.real_bodies_);

		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(water_block_inner, true);
		/**
		 * @brief 	Particle relaxation starts here.
		 */
		random_water_body_particles.exec(0.25);
		relaxation_step_inner.SurfaceBounding().exec();
		write_real_body_states.writeToFile(0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_real_body_states.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;

		/** Output results. */
		write_real_body_particle_reload_files.writeToFile(0);
		return 0;
	}
	/** external force */
	SimpleDynamics<ExternalField> drop_external_field(water_block);
	/** Initial velocity field */
	SimpleDynamics<InitialVelocity> drop_initial_velocity(water_block);
	/**
	 * @brief 	Methods used for time stepping.
	 */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannConsistency> fluid_pressure_relaxation(water_block_inner);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> fluid_density_relaxation(water_block_inner);
	InteractionWithUpdate<KernelCorrectionMatrixInner> corrected_configuration_fluid(water_block_inner, 0.5);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> fluid_density_by_summation(water_block_inner);
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, 0.0));
	SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block, gravity_ptr);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	/** We can output a method-specific particle data for debug */
	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------

	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<KineticEnergy>>>
		write_water_kinetic_energy(io_environment, water_block);
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<PotentialEnergy>>>
		write_water_potential_energy(io_environment, water_block);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);

	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	drop_initial_velocity.exec();
	//----------------------------------------------------------------------
		//	Setup for time-stepping control
		//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.RestartStep();
	int screen_output_interval = 100;
	int observation_sample_interval = 50;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 30.0; /**< End time. */
	Real D_Time = 0.01;	  /**< Time stamps for output of body states. */
	Real dt = 0.0;		  /**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	TimeInterval interval_computing_time_step;
	TimeInterval interval_computing_fluid_pressure_relaxation;
	TimeInterval interval_updating_configuration;
	TickCount time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile();
	write_water_kinetic_energy.writeToFile(number_of_iterations);
	write_water_potential_energy.writeToFile(number_of_iterations);
	write_recorded_water_pressure.writeToFile(number_of_iterations);
	/**
	 * @brief 	Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** outer loop for dual-time criteria time-stepping. */
			time_instance = TickCount::now();
			fluid_step_initialization.exec();
			Real advection_dt = 0.3 * fluid_advection_time_step.exec();
			fluid_density_by_summation.exec();
			corrected_configuration_fluid.exec();
			interval_computing_time_step += TickCount::now() - time_instance;

			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			Real acoustic_dt = 0.0;
			//while (relaxation_time < advection_dt)
			//{
				/** inner loop for dual-time criteria time-stepping.  */
				acoustic_dt = fluid_acoustic_time_step.exec();
				drop_external_field.exec();
				fluid_pressure_relaxation.exec(acoustic_dt);
				fluid_density_relaxation.exec(acoustic_dt);
				relaxation_time += acoustic_dt;
				integration_time += acoustic_dt;
				GlobalStaticVariables::physical_time_ += acoustic_dt;
			//}
			interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";
			}

			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = TickCount::now();
			water_block.updateCellLinkedList();
			water_block_inner.updateConfiguration();
			interval_updating_configuration += TickCount::now() - time_instance;
		}

		body_states_recording.writeToFile();
		write_water_kinetic_energy.writeToFile(number_of_iterations);
		write_water_potential_energy.writeToFile(number_of_iterations);
		write_recorded_water_pressure.writeToFile(number_of_iterations);
		TickCount t2 = TickCount::now();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;

	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
		<< interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
		<< interval_computing_fluid_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
		<< interval_updating_configuration.seconds() << "\n";

	return 0;
}
