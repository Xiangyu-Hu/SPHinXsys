/**
 * @file 	Turek_Hron.cpp
 * @brief 	This is the benchmark test of fluid-structure interaction.
 * @details We consider a flow-induced vibration of an elastic plate behind a cylinder in 3D.
 *			The case can be found in Chi Zhang, Massoud Rezavand, Xiangyu Hu,
 *			Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 *			Journal of Computation Physics 404 (2020) 109135.
 * @author 	Anastazja Broniatowska
 */
#include "sphinxsys.h"

#include "Turek_Hron.h" //	case file to setup the test case
using namespace SPH;

int main(int ac, char *av[])
{
	//GEOMETRY
	auto fluid_shape = createFluid();
	auto wall_shape = createWall();
	auto body_shape = createInsertBody();
	
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(wall_shape->getBounds(), resolution_ref);
	sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
	sph_system.setReloadParticles(false);		// Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
	IOEnvironment io_environment(sph_system);
	
	
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, fluid_shape);
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, wall_shape);
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();

	SolidBody insert_body(sph_system, body_shape);
	insert_body.defineAdaptationRatios(1.15, 2.0);
	insert_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	insert_body.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? insert_body.generateParticles<ParticleGeneratorReload>(io_environment, insert_body.getName())
		: insert_body.generateParticles<ParticleGeneratorLattice>();

	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_block);
	InnerRelation insert_body_inner(insert_body);
	ContactRelation water_block_contact(water_block, RealBodyVector{&wall_boundary, &insert_body});
	ComplexRelation water_block_complex(water_block_inner, water_block_contact);
	ContactRelation insert_body_contact(insert_body, {&water_block});

	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(insert_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_insert_body_to_vtp(io_environment, {&insert_body});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, {&insert_body});
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(insert_body_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_insert_body_particles.parallel_exec(0.25);
		relaxation_step_inner.SurfaceBounding().parallel_exec();
		write_insert_body_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_insert_body_to_vtp.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}

	/** Initialize particle acceleration. */
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
		// INLET
	// emmitter to inject particles
	Vec3d emitter_halfsize(0.5 * BW, 0.5 * DH, 0.5 * DH);
	Vec3d emitter_translation(0.5 * BW, 0.5 * DH, 0.5 * DH);
	BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(Transform3d(Vec3d(emitter_translation)), emitter_halfsize));
	SimpleDynamics<fluid_dynamics::EmitterInflowInjection, BodyAlignedBoxByParticle> emitter_inflow_injection(emitter, 10, 0);
	
	// inflow region to impose velocity
	Vec3d inflow_region_halfsize(0.5 * inflow_length, 0.5 * DH, 0.5 * DH);
	Vec3d inflow_region_translation(0.5 * inflow_length, 0.5 * DH, 0.5 * DH);
	BodyAlignedBoxByCell inflow_region(water_block, makeShared<AlignedBoxShape>(Transform3d(Vec3d(inflow_region_translation)), inflow_region_halfsize));
	SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>, BodyAlignedBoxByCell> emitter_buffer_inflow_condition(inflow_region);
	
	//OUTLET
	Vec3d disposer_halfsize(0.5 * BW, 0.75 * DH, 0.75 * DH);
	Vec3d disposer_translation(DL - 0.5 * BW, 0.5 * DH, 0.5 * DH);
	BodyAlignedBoxByCell disposer(water_block, makeShared<AlignedBoxShape>(Transform3d(Vec3d(disposer_translation)), disposer_halfsize));
	SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion, BodyAlignedBoxByCell> disposer_outflow_deletion(disposer, 0);
	

	/** time-space method to detect surface particles. */
	InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex>
		inlet_outlet_surface_particle_indicator(water_block_complex);

	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_complex);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	/** modify the velocity of boundary particles with free-stream velocity. */
	SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
	/** correct the velocity of boundary particles with free-stream velocity through the post process of pressure relaxation. */
	pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
	/* Density relaxation*/
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	InteractionDynamics<CombinedLocalInteraction<
		fluid_dynamics::ViscousAccelerationWithWall,
		fluid_dynamics::TransportVelocityCorrectionComplex>>
		viscous_acceleration_and_transport_correction(water_block_complex);
	/** Computing vorticity in the flow. */
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_complex.getInnerRelation());

	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(insert_body);
	/** Corrected configuration for the elastic insert body. */
	InteractionDynamics<solid_dynamics::CorrectConfiguration> insert_body_corrected_configuration(insert_body_inner);
	/** Compute the force exerted on solid body due to fluid pressure and viscosity. */
	InteractionDynamics<solid_dynamics::FluidViscousForceOnSolid> viscous_force_on_solid(insert_body_contact);
	InteractionDynamics<solid_dynamics::FluidForceOnSolidUpdate>
		fluid_force_on_solid_update(insert_body_contact, viscous_force_on_solid);
	/** Compute the average velocity of the insert body. */
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(insert_body);
	//----------------------------------------------------------------------
	//	Algorithms of solid dynamics.
	//----------------------------------------------------------------------
	/** Compute time step size of elastic solid. */
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> insert_body_computing_time_step_size(insert_body);
	/** Stress relaxation for the inserted body. */
	Dynamics1Level<solid_dynamics::Integration1stHalf> insert_body_stress_relaxation_first_half(insert_body_inner);
	Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);
	/** Constrain region of the inserted body. */
	auto constraint = createBeamBaseShape();
	BodyRegionByParticle beam_base(insert_body, constraint);
	SimpleDynamics<solid_dynamics::FixConstraint, BodyRegionByParticle> constraint_beam_base(beam_base);
	/** Update norm .*/
	SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(insert_body);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
	RegressionTestTimeAveraged<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalViscousForceOnSolid>>>
		write_total_viscous_force_on_insert_body(io_environment, insert_body);

	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	/** initialize cell linked lists for all bodies. */
	sph_system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	sph_system.initializeSystemConfigurations();
	/** computing surface normal direction for the wall. */
	wall_boundary_normal_direction.parallel_exec();
	/** computing surface normal direction for the insert body. */
	insert_body_normal_direction.parallel_exec();
	/** computing linear reproducing configuration for the insert body. */
	insert_body_corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real end_time = 200.0;
	Real output_interval = end_time / 2000.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_real_body_states.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			inlet_outlet_surface_particle_indicator.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec();

			/** FSI for viscous force. */
			viscous_force_on_solid.parallel_exec();
			/** Update normal direction on elastic body.*/
			insert_body_update_normal.parallel_exec();
			size_t inner_ite_dt = 0;
			size_t inner_ite_dt_s = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				/** Fluid pressure relaxation */
				pressure_relaxation.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_force_on_solid_update.parallel_exec();
				/** Fluid density relaxation */
				density_relaxation.parallel_exec(dt);

				/** Solid dynamics. */
				inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					Real dt_s = SMIN(insert_body_computing_time_step_size.parallel_exec(), dt - dt_s_sum);
					insert_body_stress_relaxation_first_half.parallel_exec(dt_s);
					constraint_beam_base.parallel_exec();
					insert_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
					inner_ite_dt_s++;
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				emitter_buffer_inflow_condition.parallel_exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			emitter_inflow_injection.parallel_exec();
			disposer_outflow_deletion.parallel_exec();

			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_complex.updateConfiguration();
			/** one need update configuration after periodic condition. */
			insert_body.updateCellLinkedList();
			insert_body_contact.updateConfiguration();
			/** write run-time observation into file */
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.writeToFile();
		write_total_viscous_force_on_insert_body.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	if (sph_system.generate_regression_data_)
	{
		// The lift force at the cylinder is very small and not important in this case.
		write_total_viscous_force_on_insert_body.generateDataBase({1.0e-2, 1.0e-2, 1.0e-2}, {1.0e-2, 1.0e-2, 1.0e-2});
	}
	else
	{
		write_total_viscous_force_on_insert_body.newResultTest();
	}

	return 0;
}
