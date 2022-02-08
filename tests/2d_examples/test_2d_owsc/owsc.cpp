/**
* @file 	owsc.cpp
* @brief 	This is the test of wave interaction with Oscillating Wave Surge Converter (OWSC)
* @author   Chi Zhang and Xiangyu Hu
*/
/** Header files. */
#include "sphinxsys.h"

#include "owsc.h"
/** Namespace. */
using namespace SPH;
/** The main program.*/
int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);
	/** The water block, body, material and particles container. */
	WaterBlock water_block(system, "WaterBody");
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));
	/** The wall boundary, body and particles container. */
	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	/** Flap, OWSC system. Body, material and particle container. */
	Flap flap(system, "Flap");
	ElasticSolidParticles flap_particles(flap, makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson));
	/** Pressure probe on Flap. */
	ObserverBody observer(system, "FlapObserver");
	ObserverParticles 	observer_particles(observer, makeShared<FlapObserverParticleGenerator>());
	/** topology */
	BodyRelationInner water_block_inner(water_block);
	BodyRelationInner flap_inner(flap);
	ComplexBodyRelation water_block_complex(water_block_inner, { &wall_boundary, &flap });
	BodyRelationContact flap_contact(flap, { &water_block });
	BodyRelationContact observer_contact_with_water(observer, { &water_block });
	BodyRelationContact observer_contact_with_flap(observer, {&flap});
	/**
	 * Methods only used only once
	 */
	/** corrected strong configuration. */
	solid_dynamics::CorrectConfiguration flap_corrected_configuration(flap_inner);
	/** 
	 * Methods used for time stepping
	 */
	/** Define external force.*/
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** Time step initialization, add gravity. */
	TimeStepInitialization 	initialize_gravity_to_fluid(water_block, gravity);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex);
	/** time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize	get_fluid_advection_time_step_size(water_block, U_f);
	/** time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationRiemannWithWall	pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall	density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Inflow boundary condition. */
	MultiPolygonShape damping_buffer_shape(createDampingBufferShape());
	BodyRegionByCell damping_buffer(water_block, "DampingBuffer", damping_buffer_shape);
	fluid_dynamics::DampingBoundaryCondition	damping_wave(water_block, damping_buffer);
	/** Fluid force on flap. */
	solid_dynamics::FluidForceOnSolidUpdate 	fluid_force_on_flap(flap_contact);
	/** average velocity for flap. */
	solid_dynamics::AverageVelocityAndAcceleration	average_velocity_and_acceleration(flap);
	solid_dynamics::UpdateElasticNormalDirection 	flap_update_normal(flap);
	/** constrain region of the part of wall boundary. */
	MultiPolygonShape wave_maker_shape(createWaveMakerShape());
	BodyRegionByParticle wave_maker(wall_boundary, "WaveMaker", wave_maker_shape);
	WaveMaking wave_making(wall_boundary, wave_maker);

	/** Multi-body system. */
	/** set up the multi body system. */
	SimTK::MultibodySystem MBsystem;
	/** the bodies or matter of the system. */
	SimTK::SimbodyMatterSubsystem 	matter(MBsystem);
	/** the forces of the system. */
	SimTK::GeneralForceSubsystem 	forces(MBsystem);
	/** mass properties of the fixed spot. */
	MultiPolygonShape flap_multibody_shape(createFlapSimbodyConstrainShape());
	FlapSystemForSimbody  	flap_multibody(flap, "FlapMultiBody", flap_multibody_shape);
	/** Mass properties of the constrained spot. 
	 * SimTK::MassProperties(mass, center of mass, inertia)
	 */
	SimTK::Body::Rigid    pin_spot_info(*flap_multibody.body_part_mass_properties_);
	/** 
	 * @brief   Pin (MobilizedBody &parent, const Transform &X_PF, const Body &bodyInfo, const 
	 					Transform &X_BM, Direction=Forward)1
	 * @details Create a Pin mobilizer between an existing parent (inboard) body P and 
	 * 			a new child (outboard) body B created by copying the given bodyInfo into 
	 *			a privately-owned Body within the constructed MobilizedBody object.
	 * @param[in] inboard(SimTK::Vec3) Defines the location of the joint point relative to the parent body.
	 * @param[in] outboard(SimTK::Vec3) Defines the body's origin location to the joint point. 
	 * @note	The body's origin location can be the mass center, the the center of mass should be Vec3(0)
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
     *				\a zeroHeight == 0, i.e., a body's potential energy is defined to be zero
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
	 *		also signed. 0.475
	 */
	SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, Real(-9.81), 0.0), 0.0);
	/** discrete forces acting on the bodies. */
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
	SimTK::Force::MobilityLinearDamper linear_damper(forces, pin_spot, SimTK::MobilizerUIndex(0), 20.0);
	/** Time stepping method for multibody system.*/
	SimTK::State state = MBsystem.realizeTopology();
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	integ.initialize(state);
	/**
	* Coupling between SimBody and SPH.
	*/
	solid_dynamics::TotalForceOnSolidBodyPartForSimBody
		force_on_spot_flap(flap, flap_multibody, MBsystem, pin_spot, force_on_bodies, integ);
	solid_dynamics::ConstrainSolidBodyPartBySimBody
		constraint_spot_flap(flap, flap_multibody, MBsystem, pin_spot, force_on_bodies, integ);
	/** Output. */
	In_Output in_output(system);
	BodyStatesRecordingToVtp 		write_real_body_states(in_output, system.real_bodies_);
	RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<solid_dynamics::TotalForceOnSolid>> write_total_force_on_flap(in_output, flap);
	WriteSimBodyPinData			write_flap_pin_data(in_output, integ, pin_spot);
	/** WaveProbes. */
	MultiPolygonShape wave_probe_buffer_no_4_shape(createWaveProbeShape4());
	BodyRegionByCell wave_probe_buffer_no_4(water_block, "WaveProbe_04", wave_probe_buffer_no_4_shape);
	BodyReducedQuantityRecording<fluid_dynamics::FreeSurfaceProbeOnFluidBody>
		wave_probe_4(in_output, water_block,  wave_probe_buffer_no_4);
	
	MultiPolygonShape wave_probe_buffer_no_5_shape(createWaveProbeShape5());
	BodyRegionByCell wave_probe_buffer_no_5(water_block, "WaveProbe_05", wave_probe_buffer_no_5_shape);
	BodyReducedQuantityRecording<fluid_dynamics::FreeSurfaceProbeOnFluidBody>
		wave_probe_5(in_output, water_block,  wave_probe_buffer_no_5);
	
	MultiPolygonShape wave_probe_buffer_no_12_shape(createWaveProbeShape12());
	BodyRegionByCell wave_probe_buffer_no_12(water_block, "WaveProbe_12", wave_probe_buffer_no_12_shape);
	BodyReducedQuantityRecording<fluid_dynamics::FreeSurfaceProbeOnFluidBody>
		wave_probe_12(in_output, water_block, wave_probe_buffer_no_12);
	/** Pressure probe. */
	ObservedQuantityRecording<Real> pressure_probe("Pressure", in_output, observer_contact_with_water);
	/** Interpolate the particle position in flap to move the observer accordingly. */
	observer_dynamics::InterpolatingAQuantity<Vecd>
		interpolation_observer_position(observer_contact_with_flap, "Position", "Position");
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	/** offset particle position */
	flap_particles.offsetInitialParticlePosition(offset);
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	flap_particles.initializeNormalDirectionFromBodyShape();
	flap_corrected_configuration.parallel_exec();

	write_real_body_states.writeToFile(0);
	write_total_force_on_flap.writeToFile(0);
	write_flap_pin_data.writeToFile(0);
	wave_probe_4.writeToFile(0);
	wave_probe_5.writeToFile(0);
	wave_probe_12.writeToFile(0);
	pressure_probe.writeToFile(0);
	/** Simulation start here. */
	/** starting time zero. */
	GlobalStaticVariables::physical_time_ = 0.0;
	int number_of_iterations = 0;
	int screen_output_interval = 1000;
	Real End_Time = total_physical_time;
	Real D_Time = End_Time/100.0;
	Real dt = 0.0; 
	Real total_time = 0.0;
	Real relax_time = 1.0;
	/** statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main Loop. */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integral_time = 0.0;
		while (integral_time < D_Time) 
		{
			initialize_gravity_to_fluid.parallel_exec();
			
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			/** Viscous force exerting on flap. */
			fluid_force_on_flap.viscous_force_.parallel_exec();
			flap_update_normal.parallel_exec();
			
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{
				pressure_relaxation.parallel_exec(dt);
				fluid_force_on_flap.parallel_exec();
				density_relaxation.parallel_exec(dt);
				/** solid dynamics. */
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				if(total_time >= relax_time)
				{
					SimTK::State& state_for_update = integ.updAdvancedState();
					Real angle = pin_spot.getAngle(state_for_update);
					force_on_bodies.clearAllBodyForces(state_for_update);
					force_on_bodies.setOneBodyForce(state_for_update, pin_spot, force_on_spot_flap.parallel_exec(angle));
					integ.stepBy(dt);
					constraint_spot_flap.parallel_exec();
					wave_making.parallel_exec(dt);
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);
				interpolation_observer_position.parallel_exec();

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integral_time += dt;
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
			damping_wave.parallel_exec(Dt);
			water_block.updateCellLinkedList();
			wall_boundary.updateCellLinkedList();
			flap.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			flap_contact.updateConfiguration();
			observer_contact_with_water.updateConfiguration();
			if(total_time >= relax_time)
			{
				write_total_force_on_flap.writeToFile(number_of_iterations);
				write_flap_pin_data.writeToFile(GlobalStaticVariables::physical_time_ );
				wave_probe_4.writeToFile(number_of_iterations);
				wave_probe_5.writeToFile(number_of_iterations);
				wave_probe_12.writeToFile(number_of_iterations);
				pressure_probe.writeToFile(number_of_iterations);
			}
		}

		tick_count t2 = tick_count::now();
		if(total_time >= relax_time)
			write_real_body_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	
	write_total_force_on_flap.newResultTest();

	return 0;
}
