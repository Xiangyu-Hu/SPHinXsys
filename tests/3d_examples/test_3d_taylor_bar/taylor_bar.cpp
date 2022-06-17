/**
 * @file taylor_bar.cpp
 * @brief This is the case setup for plastic taylor bar.
 * @author Chi Zhang  and Xiangyu Hu
 * @ref 	doi.org/10.1007/s40571-019-00277-6
 */
#include "sphinxsys.h"

#include "taylor_bar.h" /**< Case setup for this example. */

using namespace SPH;

int main(int ac, char *av[])
{
	/** Setup the system. Please the make sure the global domain bounds are correctly defined. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);
	system.run_particle_relaxation_ = false;
	system.reload_particles_ = true;
// handle command line arguments
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	InOutput in_output(system);

	/** create a body with corresponding material, particles and reaction model. */
	SolidBody column(system, makeShared<Column>("Column"));
	column.sph_adaptation_->resetAdaptationRatios(1.3, 1.0);
	column.defineBodyLevelSetShape()->writeLevelSet(column);
	column.defineParticlesAndMaterial<ElasticSolidParticles, HardeningPlasticSolid>(
		rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);
	(!system.run_particle_relaxation_ && system.reload_particles_)
		? column.generateParticles<ParticleGeneratorReload>(in_output, column.getBodyName())
		: column.generateParticles<ParticleGeneratorLattice>();
	column.addBodyStateForRecording<Vecd>("NormalDirection");

	SolidBody wall(system, makeShared<Wall>("Wall"));
	wall.defineParticlesAndMaterial<SolidParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
	wall.generateParticles<ParticleGeneratorLattice>();

	/** Define Observer. */
	ObserverBody my_observer(system, "MyObserver");
	my_observer.generateParticles<ColumnObserverParticleGenerator>();

	/**body relation topology */
	BodyRelationInner column_inner(column);
	BodyRelationContact my_observer_contact(my_observer, {&column});
	SolidBodyRelationContact column_wall_contact(column, {&wall});
	/**define simple data file input and outputs functions. */
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);

	if (system.run_particle_relaxation_)
	{
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		/** Random reset the insert body particle position. */
		RandomizeParticlePosition random_column_particles(column);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_column_to_vtp(in_output, column);
		/** Write the particle reload files. */

		ReloadParticleIO write_particle_reload_files(in_output, {&column});
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(column_inner);
		/**
		 * @brief 	Particle relaxation starts here.
		 */
		random_column_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_states.writeToFile(0.0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << setprecision(9) << "Relaxation steps for the column body N = " << ite_p << "\n";
				write_column_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of cylinder body finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.writeToFile(0.0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	All numerical methods will be used in this case.
	//----------------------------------------------------------------------
	SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall);
	InitialCondition initial_condition(column);
	/** Corrected configuration. */
	solid_dynamics::CorrectConfiguration corrected_configuration(column_inner);
	/** Time step size calculation with given CFL number. */
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(column, 0.3);

	/** stress and deformation relaxation. */
	solid_dynamics::PlasticStressRelaxationFirstHalf stress_relaxation_first_half(column_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(column_inner);
	solid_dynamics::DynamicContactForceWithWall column_wall_contact_force(column_wall_contact);

	//----------------------------------------------------------------------
	//	Output
	//----------------------------------------------------------------------
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_velocity("Velocity", in_output, my_observer_contact);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_displacement("Position", in_output, my_observer_contact);

	//----------------------------------------------------------------------
	// From here the time stepping begins.
	//----------------------------------------------------------------------
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_normal_direction.parallel_exec();
	corrected_configuration.parallel_exec();
	initial_condition.parallel_exec();
	//----------------------------------------------------------------------
	// Setup time-stepping related simulation parameters.
	//----------------------------------------------------------------------
	int ite = 0;
	Real end_time = 1.0e-4;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	Real output_period = 1.0e-6; // anyway 50 write_states files in total
	Real dt = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_states.writeToFile();
	write_displacement.writeToFile(0);
	write_velocity.writeToFile(0);
	//----------------------------------------------------------------------
	// Main time-stepping loop.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			if (ite % screen_output_interval == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";

				if (ite != 0 && ite % observation_sample_interval == 0)
				{
					write_displacement.writeToFile(ite);
					write_velocity.writeToFile(ite);
				}
			}
			column_wall_contact_force.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			column.updateCellLinkedList();
			column_wall_contact.updateConfiguration();

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	write_displacement.newResultTest();
	write_velocity.newResultTest();

	return 0;
}
