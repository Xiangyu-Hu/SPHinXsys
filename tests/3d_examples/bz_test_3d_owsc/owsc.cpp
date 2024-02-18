/**
* @file 	owsc.cpp
* @brief 	This is the test of wave interaction with Oscillating Wave Surge Converter (OWSC) in 3D wave tank.
* @author   Chi Zhang and Xiangyu Hu
* @version  0.3.0
* @note  	Observer, moving with mobile solid body, can find template in this case.
*			-- Chi ZHANG
*/
/** header file and namespace. */
#include "sphinxsys.h"
/** case file to setup the test case */
#include "case.h"
/** Name space. */
using namespace SPH;
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, dp_0);
	IOEnvironment io_environment(system);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.setRunParticleRelaxation(false);
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.setReloadParticles(true);

	/** The water block, body, material and particles container. */
	FluidBody water_block(system, makeShared<Water>("WaterBody"));
	water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	water_block.addBodyStateForRecording<Real>("Pressure");

	/** wall. */
	SolidBody wall_boundary(system, makeShared<MyTank>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();

	SolidBody wave_maker(system, makeShared<MyWaveMaker>("WaveMaker"));;
	wave_maker.defineParticlesAndMaterial<SolidParticles, Solid>();
	wave_maker.generateParticles<ParticleGeneratorLattice>();

	/** The flap. */
	SolidBody flap(system, makeShared<MyFlap>("Flap"));
	flap.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s);
	flap.generateParticles<ParticleGeneratorLattice>();

	/** Pressure probe on Flap. */
	ObserverBody observer(system, "FlapObserver");
	observer.generateParticles<FlapObserverParticleGenerator>();

	/** Gravity. */
	Gravity gravity(Vec3d(0.0,-gravity_g, 0.0));
	/** topology */
	InnerRelation water_block_inner(water_block);
	InnerRelation flap_inner(flap);
	ComplexRelation water_block_complex(water_block_inner, { &wall_boundary, &wave_maker,&flap });
	ComplexRelation wave_maker_complex(wave_maker, { &water_block });
	ComplexRelation flap_complex(flap, { &water_block });
	ContactRelation flap_contact(flap, { &water_block });
	ContactRelation observer_contact_with_water(observer, { &water_block });
	ContactRelation observer_contact_with_flap(observer, { &flap });

	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.RunParticleRelaxation())
	{
		InnerRelation wall_inner(wall_boundary);
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_flap_particles(flap);
		SimpleDynamics<RandomizeParticlePosition> random_wall_particles(wall_boundary);
		/** Write the body state to Vtu file. */
		BodyStatesRecordingToVtp write_relaxed_body(io_environment, { &flap, &wall_boundary});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, { &flap, &wall_boundary });
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_flap_inner(flap_inner);
		relax_dynamics::RelaxationStepInner relaxation_wall_inner(wall_inner);
		/** @brief 	Particle relaxation starts here.*/
		random_flap_particles.exec(0.25);
		random_wall_particles.exec(0.25);
		relaxation_flap_inner.SurfaceBounding().exec();
		relaxation_wall_inner.SurfaceBounding().exec();
		write_relaxed_body.writeToFile(0);
		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_flap_inner.exec();
			relaxation_wall_inner.exec();
			ite_p += 1;
			if (ite_p % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_relaxed_body.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	/**
	 * Methods only used only once
	 */
	SimpleDynamics<OffsetInitialPosition> flap_offset_position(flap, offset);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> wave_maker_normal_direction(wave_maker);
	SimpleDynamics<NormalDirectionFromBodyShape> flap_normal_direction(flap);

	InteractionWithUpdate<KernelCorrectionMatrixComplex> corrected_configuration_fluid(water_block_complex, 0.95);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> corrected_configuration_wave_maker(wave_maker_complex);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> corrected_configuration_flap(flap_complex);
	/** Time step initialization, add gravity. */
	SimpleDynamics<TimeStepInitialization> initialize_time_step_to_fluid(water_block, makeShared<Gravity>(Vec3d(0.0, -gravity_g, 0.0)));
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_density_by_summation(water_block_complex);
	/** time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	/** time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	/** pressure relaxation using verlet time stepping. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannConsistencyWithWall> pressure_relaxation(water_block_complex);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
	/** Inflow boundary condition. */
	//fluid_dynamics::DampingBoundaryCondition	damping_wave(my_water, new DampingBuffer(my_water, "DampingBuffer"), 5.0);
	BodyRegionByCell damping_buffer(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_damping), halfsize_damping));
	SimpleDynamics<fluid_dynamics::DampingBoundaryCondition> damping_wave(damping_buffer);
	/** Update pressure for fluid particles after sorting. */
	//fluid_dynamics::UpdatePressureOfSortingParticleForObserving update_pressure(my_water);
	/** Fluid froce on flap. */
	InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(flap_contact);
	InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_flap(flap_contact, viscous_force_on_solid);
	/** constrain region of the part of wall boundary. */
	SimpleDynamics<WaveMaking> wave_making(wave_maker);

	/** Muti-body system. */
	/** set up the multi body system. */
	SimTK::MultibodySystem MBsystem;
	/** the bodies or matter of the system. */
	SimTK::SimbodyMatterSubsystem 	matter(MBsystem);
	/** the forces of the system. */
	SimTK::GeneralForceSubsystem 	forces(MBsystem);
	/** mass proeprties of the fixed spot. */
	//FlapSystemForSimbody*  	flap_multibody = new FlapSystemForSimbody(my_flap, "OWSC-Pin", rho0_s);
	FlapSystemForSimbody flap_multibody(flap, makeShared<TriangleMeshShapeSTL>(path_to_flap_stl, translation_flap, length_scale, "FlapMultiBody"));
	//FlapSystemForSimbody flap_multibody(flap, makeShared<TransformShape<GeometricShapeBox>>(CreateCADGeometryForOWSC(), "FlapMultiBody"));
	/** Mass properties of the consrained spot. 
	 * SimTK::MassProperties(mass, center of mass, inertia)
	 */
	SimTK::Body::Rigid   pin_spot_info(*flap_multibody.body_part_mass_properties_);
	/** 
	 * @brief   Pin (MobilizedBody &parent, const Transform &X_PF, const Body &bodyInfo, const 
	 					Transform &X_BM, Direction=Forward)1
	 * @details Create a Pin mobilizer between an existing parent (inboard) body P and 
	 * 			a new child (outboard) body B created by copying the given bodyInfo into 
	 *			a privately-owned Body within the constructed MobilizedBody object.
	 * @param[in] inboard(SimTK::Vec3) Defines the location of the joint point relative to the parent body.
	 * @param[in] outboard(SimTK::Vec3) Defines the body's origin location to the jiont point. 
	 * @note	The body's rogin location can be the mass center, the the center of mass shouled be Vec3(0)
	 * 			in SimTK::MassProperties(mass, com, inertia)
	 */
	SimTK::MobilizedBody::Pin pin_spot(matter.Ground(), SimTK::Transform(SimTK::Vec3(7.92, 0.315, 0.0)), 
		pin_spot_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
	/** set the default angle of the pin. */
	pin_spot.setDefaultAngle(0);
	/** 
	 * @details Add gravity to mb body. 
	 * @param[in,out] forces, The subsystem to which this force should be added.
	 * @param[in]     matter, The subsystem containing the bodies that will be affected.
	 * @param[in]    gravity, The default gravity vector v, interpreted as v=g*d where g=|\a gravity| is 
     *				a positive scalar and d is the "down" direction unit vector d=\a gravity/g.
	 * @param[in]  zeroHeight This is an optional specification of the default value for the height
     *				up the gravity vector that is considered to be "zero" for purposes of
     *				calculating the gravitational potential energy. The default is 
     *				 zeroHeight == 0, i.e., a body's potential energy is defined to be zero
     *				when the height of its mass center is the same as the height of the Ground
     *				origin. The zero height will have the value specified here unless
     *				explicitly changed within a particular State use the setZeroHeight() 
     *			method. 
	 * @par Force Each body B that has not been explicitly excluded will experience a force 
	 *		fb = mb*g*d, applied to its center of mass, where mb is the mass of body B. 
	 * @par Potential Energy
	 *		Gravitational potential energy for a body B is mb*g*hb where hb is the height of 
	 *		body B's mass center over an arbitrary "zero" height hz (default is hz=0), 
	 *		measured along the "up" direction -d. If pb is the Ground frame vector giving 
	 *		the position of body B's mass center, its height over or under hz is 
	 *		hb=pb*(-d) - hz. Note that this is a signed quantity so the potential energy is 
	 *		also signed. 
	 */
	SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, -gravity_g, 0.0), 0.0);
	/** discreted forces acting on the bodies. */
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	/**
	 * Add a linear damping force to the mobilized body.
	 * @class SimTK::Force::MobilityLinearDamper::MobilityLinearDamper( 	
	 * @param[in]	GeneralForceSubsystem &  	forces,
	 * @param[in]	const MobilizedBody &  	mobod,
     * @param[in]	MobilizerUIndex  whichU, e.g., MobilizerUIndex(0)
	 * @param[in]	Real  	Dampingconstant ) 
	 * Here, The damping constant c is provided, with the generated force being -c*u where u is the mobility's generalized speed.	
	 */
	SimTK::Force::MobilityLinearDamper linear_damper(forces, pin_spot, SimTK::MobilizerUIndex(0), 0.0);
	/** Time steping method for multibody system.*/
	SimTK::State state = MBsystem.realizeTopology();
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	integ.initialize(state);
	/**
	* Coupling between SimBody and SPH.
	*/
	ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
		force_on_spot_flap(flap_multibody, MBsystem, pin_spot, integ);
	SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
		constraint_spot_flap(flap_multibody, MBsystem, pin_spot, integ);

	BodyStatesRecordingToPlt write_real_body_states(io_environment, system.real_bodies_);
	RegressionTestDynamicTimeWarping<
		ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>> write_total_force_on_flap(io_environment, fluid_force_on_flap, "TotalForceOnSolid");
	WriteSimBodyPinData write_flap_pin_data(io_environment, integ, pin_spot);

	/** WaveProbes. */
	BodyRegionByCell wave_probe_buffer_no_4(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_No4), halfsize_No4, "WaveProbe_04"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
		wave_probe_4(io_environment, wave_probe_buffer_no_4);

	BodyRegionByCell wave_probe_buffer_no_5(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_No5), halfsize_No5,  "WaveProbe_05"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
		wave_probe_5(io_environment, wave_probe_buffer_no_5);

	BodyRegionByCell wave_probe_buffer_no_12(water_block, makeShared<TransformShape<GeometricShapeBox>>(Transform(translation_No12), halfsize_No12, "WaveProbe_12"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
		wave_probe_12(io_environment, wave_probe_buffer_no_12);

	/** Pressure probe. */
	ObservedQuantityRecording<Real> pressure_probe("Pressure", io_environment, observer_contact_with_water);
	// Interpolate the particle position in flap to move the observer accordingly. 
	// Seems not used? TODO: observe displacement more accurate.
	InteractionDynamics<InterpolatingAQuantity<Vecd>>
		interpolation_observer_position(observer_contact_with_flap, "Position", "Position");
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	flap_offset_position.exec();
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_normal_direction.exec();
	wave_maker_normal_direction.exec();
	flap_normal_direction.exec();

	write_real_body_states.writeToFile(0);
	write_total_force_on_flap.writeToFile(0);
	write_flap_pin_data.writeToFile(0);
	wave_probe_4.writeToFile(0);
	wave_probe_5.writeToFile(0);
	wave_probe_12.writeToFile(0);
	pressure_probe.writeToFile(0);
	/** starting simulation */
	int number_of_iterations = 0;
	int screen_output_interval = 100;
	Real relax_time = 1.0;
	Real End_Time = total_physical_time + relax_time;
	Real D_Time = total_physical_time / 500.0;
	Real Observe_Time = 0.01;
	Real Dt = 0.0;
	Real dt = 0.0; 
	Real total_time = 0.0;
	/** statistics for computing time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	/** Main Loop. */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < D_Time) 
		{
			Real observ_out_time = 0.0;
			while (observ_out_time < Observe_Time) 
			{ 
				initialize_time_step_to_fluid.exec();
				Dt = get_fluid_advection_time_step_size.exec();
				update_density_by_summation.exec();
				viscous_acceleration.exec();
				/** Viscous force exerting on flap. */
				viscous_force_on_solid.exec();
				corrected_configuration_fluid.exec();
				corrected_configuration_wave_maker.exec();
				corrected_configuration_flap.exec();

				Real relaxation_time = 0.0;
				while (relaxation_time < Dt) 
				{
					pressure_relaxation.exec(dt);
					fluid_force_on_flap.exec();
					density_relaxation.exec(dt);
			
					if(total_time >= relax_time)
					{
						SimTK::State& state_for_update = integ.updAdvancedState();
						Real angle = pin_spot.getAngle(state_for_update);
						force_on_bodies.clearAllBodyForces(state_for_update);
						force_on_bodies.setOneBodyForce(state_for_update, pin_spot, force_on_spot_flap.exec(angle));
						integ.stepBy(dt);
						constraint_spot_flap.exec();
						wave_making.exec(dt);
					}
					
					interpolation_observer_position.exec();

					dt = get_fluid_time_step_size.exec();
					relaxation_time += dt;
					integeral_time += dt;
					observ_out_time += dt;
					total_time += dt;
					if(total_time >= relax_time)
						GlobalStaticVariables::physical_time_ += dt;
				}
				
				if (number_of_iterations % screen_output_interval == 0)
				{
					std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
						<< "	Total Time = " << total_time
						<< "	Physical Time = " << GlobalStaticVariables::physical_time_
						<< "	Dt = " << Dt << "	dt = " << dt << "\n";
				}
				number_of_iterations++;
				damping_wave.exec(Dt);
				water_block.updateCellLinkedListWithParticleSort(100);
				wall_boundary.updateCellLinkedList();
				wave_maker.updateCellLinkedList();
				flap.updateCellLinkedList();
				water_block_complex.updateConfiguration();
				wave_maker_complex.updateConfiguration();
				flap_complex.updateConfiguration();
				flap_contact.updateConfiguration();
				observer_contact_with_water.updateConfiguration();

			}
			if(total_time >= relax_time)
			{
				write_total_force_on_flap.writeToFile(number_of_iterations);
				write_flap_pin_data.writeToFile(GlobalStaticVariables::physical_time_);
				wave_probe_4.writeToFile(number_of_iterations);
				wave_probe_5.writeToFile(number_of_iterations);
				wave_probe_12.writeToFile(number_of_iterations);
				pressure_probe.writeToFile(number_of_iterations);
			}
		}

		TickCount t2 = TickCount::now();
		if (total_time >= relax_time)
			write_real_body_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
