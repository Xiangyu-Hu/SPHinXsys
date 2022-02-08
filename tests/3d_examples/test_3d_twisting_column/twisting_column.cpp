/**
 * @file twisting_column.cpp
 * @brief This is an example of solid with classic neohookean model 
 * to demonstrate the robustness of the formulation with Kirchhoff stress decomposition.
 * @author Chi Zhang  and Xiangyu Hu
 * @ref 	DOI: 10.1016/j.cma.2014.09.024
 */
#include "sphinxsys.h"

#include "twisting_column.h" /**< Case setup for this example. */

using namespace SPH;

int main()
{
	/** Setup the system. Please the make sure the global domain bounds are correctly defined. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);
	/** create a body with corresponding material, particles and reaction model. */
	Column column(system, "Column");
	ElasticSolidParticles particles(column, makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson));
	/** Define Observer. */
	ObserverBody my_observer(system, "MyObserver");
	ObserverParticles observer_particles(my_observer, makeShared<ObserverParticleGenerator>());
	/**body relation topology */
	BodyRelationInner column_inner(column);
	BodyRelationContact my_observer_contact(my_observer, {&column});
	//----------------------------------------------------------------------
	//	All numerical methods will be used in this case.
	//----------------------------------------------------------------------
	InitialCondition initial_condition(column);
	/** Corrected configuration. */
	solid_dynamics::CorrectConfiguration corrected_configuration(column_inner);
	/** Time step size calculation. We use CFL = 0.5 due to the very large twisting speed. */
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(column, 0.5);
	/** active and passive stress relaxation. */
	solid_dynamics::KirchhoffStressRelaxationFirstHalf stress_relaxation_first_half(column_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(column_inner);
	/** Constrain the holder. */
	TriangleMeshShapeBrick holder_shape(halfsize_holder, resolution, translation_holder);
	BodyRegionByParticle holder(column, "Holder", holder_shape);
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(column, holder);
	//----------------------------------------------------------------------
	//	Output
	//----------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_velocity("Velocity", in_output, my_observer_contact);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_displacement("Position", in_output, my_observer_contact);
	//----------------------------------------------------------------------
	// From here the time stepping begins.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	initial_condition.parallel_exec();
	corrected_configuration.parallel_exec();
	write_states.writeToFile(0);
	write_displacement.writeToFile(0);
	write_velocity.writeToFile(0);
	//----------------------------------------------------------------------
	// Setup time-stepping related simulation parameters.
	//----------------------------------------------------------------------
	int ite = 0;
	Real end_time = 0.5;
	Real output_period = end_time / 250.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	// Main time-stepping loop.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
			write_displacement.writeToFile(ite);
			write_velocity.writeToFile(ite);
		}
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_displacement.newResultTest();
	write_velocity.newResultTest();

	return 0;
}
