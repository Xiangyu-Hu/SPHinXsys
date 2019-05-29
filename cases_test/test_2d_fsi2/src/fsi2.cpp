/**
 * @file 	fsi2.cpp
 * @brief 	This is the benchmark test of fliud-structure interaction.
 * @details We consider a flow-induced vibration of an elastic beam behind a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 * @version 0.1
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
Real DL = 11.0; 					/**< Channel length. */
Real DH = 4.1; 						/**< Channel height. */
Real particle_spacing_ref = 0.1; 	/**< Initial reference particle spacing. */

Real DLsponge = particle_spacing_ref *20.0;	/**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0; 		/**< Boundary width, determined by spcific layer of boundary patciels. */

Vec2d insert_circle_center(2.0, 2.0);		/**< Location of the cylinder center. */
Real insert_circle_radius = 0.5;			/**< Radius of the cylinder. */

Real bh = 0.4*insert_circle_radius;			/**< Height of the beam. */
Real bl = 7.0*insert_circle_radius;			/**< Length of the beam. */
/**
 * @brief Geomerty of the beam. Defined through the 4 corners of a box.
 */
Real hbh = bh / 2.0;
Vec2d BLB(insert_circle_center[0], insert_circle_center[1] - hbh);
Vec2d BLT(insert_circle_center[0], insert_circle_center[1] + hbh);
Vec2d BRB(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] - hbh);
Vec2d BRT(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] + hbh);
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		/**< Density. */
Real U_f = 1.0;			/**< Cheractristic velocity. */
Real c_f = 10.0*U_f;	/**< Speed of sound. */
Real Re = 100.0;		/**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re;	/**< Dynamics visocisty. */
Real k_f = 0.0;			/**< kinetic smoothness. */

Real initial_pressure = 0.0;		/**< Initial pressure of fluid. */
Vec2d intial_velocity(0.0, 0.0);	/**< Initial velocity of fluid. */
/**
 * @brief Material properties of the solid,
 */
Real rho0_s = 10.0; 		/**< Reference density.*/
Real poisson = 0.4; 		/**< Poisson ratio.*/
Real Ae = 1.4e3; 			/**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
/**
* @brief define geometry and initial conditions of SPH bodies
*/
/**
* @brief create a water block shape
*/
std::vector<Point> CreatWaterBlockShape()
{
	//geometry
	std::vector<Point> water_block_shape;
	water_block_shape.push_back(Point(-DLsponge, 0.0));
	water_block_shape.push_back(Point(-DLsponge, DH));
	water_block_shape.push_back(Point(DL, DH));
	water_block_shape.push_back(Point(DL, 0.0));
	water_block_shape.push_back(Point(-DLsponge, 0.0));

	return water_block_shape;
}
/**
* @brief create a water block buffer shape
*/
std::vector<Point> CreatInflowBufferShape()
{
	std::vector<Point> inlfow_buffer_shape;
	inlfow_buffer_shape.push_back(Point(-DLsponge, 0.0));
	inlfow_buffer_shape.push_back(Point(-DLsponge, DH));
	inlfow_buffer_shape.push_back(Point(0.0, DH));
	inlfow_buffer_shape.push_back(Point(0.0, 0.0));
	inlfow_buffer_shape.push_back(Point(-DLsponge, 0.0));

	return inlfow_buffer_shape;
}
/**
* @brief create a beam shape
*/
std::vector<Point> CreatBeamShape()
{
	std::vector<Point> beam_shape;
	beam_shape.push_back(BLB);
	beam_shape.push_back(BLT);
	beam_shape.push_back(BRT);
	beam_shape.push_back(BRB);
	beam_shape.push_back(BLB);

	return beam_shape;
}
/**
* @brief create outer wall shape
*/
std::vector<Point> CreatOuterWallShape()
{
	std::vector<Point> outer_wall_shape;
	outer_wall_shape.push_back(Point(-DLsponge - BW, -BW));
	outer_wall_shape.push_back(Point(-DLsponge - BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, -BW));
	outer_wall_shape.push_back(Point(-DLsponge - BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Point> CreatInnerWallShape()
{
	std::vector<Point> inner_wall_shape;
	inner_wall_shape.push_back(Point(-DLsponge - 2.0*BW, 0.0));
	inner_wall_shape.push_back(Point(-DLsponge - 2.0*BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0*BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0*BW, 0.0));
	inner_wall_shape.push_back(Point(-DLsponge - 2.0*BW, 0.0));

	return inner_wall_shape;
}
/**
 * @brief Fluid body definition. 
 */
class WaterBlock : public WeaklyCompressibleFluidBody
{
	public:
		WaterBlock(SPHSystem &system, string body_name,
			WeaklyCompressibleFluid* material, 
			WeaklyCompressibleFluidParticles &weakly_compressible_fluid_particles, 
			int refinement_level, ParticlesGeneratorOps op)
			: WeaklyCompressibleFluidBody(system, body_name, material, 
				weakly_compressible_fluid_particles, refinement_level, op)
		{
			/** Geomerty definition. */
			body_region_.add_geometry(new Geometry(CreatWaterBlockShape()), RegionBooleanOps::add);
			/** Geomerty definition. */
			body_region_.add_geometry(new Geometry(insert_circle_center, insert_circle_radius, 100), RegionBooleanOps::sub);
			body_region_.add_geometry(new Geometry(CreatBeamShape()), RegionBooleanOps::sub);
			/** Finalize the geometry definition and correspoding opertation. */
			body_region_.done_modeling();
		}
	/**
	 * @brief Definiton of the initial condition for starting simulation. 
	 */
	void InitialCondition() 
	{
		SetAllParticleAtRest();
	}
};
/**
 * @brief Definition of the solid body.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name, 
		SolidBodyParticles &solid_particles,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, solid_particles, refinement_level, op)
	{	
		/** Geomerty definition. */
		body_region_.add_geometry(new Geometry(CreatOuterWallShape()), RegionBooleanOps::add);
		/** Geomerty definition. */
		body_region_.add_geometry(new Geometry(CreatInnerWallShape()), RegionBooleanOps::sub);
		/** Finalize the geometry definition and correspoding opertation. */
		body_region_.done_modeling();

	}
	/**
	 * @brief Definiton of the initial condition for starting simulation. 
	 */
	void InitialCondition() 
	{
		SetAllParticleAtRest();
	}
};

//insert elastic body with constraint
/**
 * @brief Definition of the inserted body as a elastic structure.
 */
class InsertedBody : public ElasticBody
{
public:
	InsertedBody(SPHSystem &system, string body_name, ElasticSolid* material,
		ElasticBodyParticles &elastic_particles,
		int refinement_level, ParticlesGeneratorOps op)
		: ElasticBody(system, body_name, 
			material, elastic_particles, refinement_level, op)
	{
		/** Geomerty definition. */
		Geometry *circle_geometry = new Geometry(insert_circle_center, insert_circle_radius, 100);
		body_region_.add_geometry(circle_geometry, RegionBooleanOps::add);
		Geometry *beam_geometry = new Geometry(CreatBeamShape());
		body_region_.add_geometry(beam_geometry, RegionBooleanOps::add);
		/** Finalize the geometry definition and correspoding opertation. */
		body_region_.done_modeling();
	}
	/**
	 * @brief Definiton of the initial condition for starting simulation. 
	 */
	void InitialCondition()
	{
		SetAllParticleAtRest();
	}
};
/**
* @brief constrain the beam base
*/
class BeamBase : public SolidBodyPart
{
public:
	BeamBase(SolidBody *solid_body, string constrianed_region_name)
		: SolidBodyPart(solid_body, constrianed_region_name)
	{
		//geometry
		Geometry *circle_geometry = new Geometry(insert_circle_center, insert_circle_radius, 100);
		soild_body_part_region_.add_geometry(circle_geometry, RegionBooleanOps::add);
		Geometry * beam_gemetry = new Geometry(CreatBeamShape());
		soild_body_part_region_.add_geometry(beam_gemetry, RegionBooleanOps::sub);
		//finish the region modeling
		soild_body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
	}
};
/**
* @brief inflow buffer
*/
class InflowBuffer : public FluidBodyPart
{
public:
	InflowBuffer(WeaklyCompressibleFluidBody* fluid_body, string constrianed_region_name)
		: FluidBodyPart(fluid_body, constrianed_region_name)
	{
		/** Geomerty definition. */
		fluid_body_part_region_.add_geometry(new Geometry(CreatInflowBufferShape()), RegionBooleanOps::add);
		/** Finalize the geometry definition and correspoding opertation. */
		fluid_body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartCells();
	}
};
/**
 * @brief Definition of an observer body with one particle located at specific position.
 */
class Observer : public ObserverLagrangianBody
{
public:
	Observer(SPHSystem &system, string body_name,
		ObserverParticles &observer_particles,
		int refinement_level, ParticlesGeneratorOps op)
		: ObserverLagrangianBody(system, body_name, observer_particles, refinement_level, op)
	{
		/** postion and volume. */
		body_input_points_volumes_.push_back(make_pair(0.5 * (BRT + BRB), 0.0));
	}
};
/**
 * @brief Inflow BCs.
 */
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(WeaklyCompressibleFluidBody* fluid_body, 
		FluidBodyPart *constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region)
	{
		u_ave_ = 0.0;
		u_ref_ = 1.0;
		t_ref = 2.0;
	}
	Vecd GetInflowVelocity(Vecd &position, Vecd &velocity) 
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0) {
			u = 6.0*u_ave_*position[1] * (DH - position[1]) / DH / DH;
			v = 0.0;
		}
		return Vecd(u, v);
	}
	void PrepareConstraint() override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5*u_ref_*(1.0 - cos(pi*run_time/ t_ref)) : u_ref_;
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
	SPHSystem system(Vec2d(-DLsponge - BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Set the starting time. */
	/**
	 * @brief Note that in this case, BCs(inflow) is time dependent, please
	 *			set physical_time_ properly when you using resart technique.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.restart_step_ = 0;
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	SymmetricTaitFluid 					fluid("Water", rho0_f, c_f, mu_f, k_f);
	WeaklyCompressibleFluidParticles 	fluid_particles("WaterBody");
	WaterBlock *water_block  = 			new WaterBlock(system, "WaterBody", &fluid, 
											fluid_particles, 0, ParticlesGeneratorOps::lattice);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	SolidBodyParticles 				solid_particles("Wall");
	WallBoundary *wall_boundary = 	new WallBoundary(system, "Wall", solid_particles,	
											0, ParticlesGeneratorOps::lattice);
	/**
	 * @brief 	Material property, particle and body creation of elastic beam(inserted body).
	 */
	ElasticSolid 					solid_material("ElasticSolid", rho0_s, Youngs_modulus, poisson, 0.0);
	ElasticBodyParticles 			inserted_body_particles("InsertedBody");
	InsertedBody *inserted_body = 	new InsertedBody(system, "InsertedBody", &solid_material, 
											inserted_body_particles, 1, ParticlesGeneratorOps::lattice);
	/**
	 * @brief 	Particle and body creation of gate observer.
	 */
	ObserverParticles 			observer_particles("Observer");
	Observer *beam_observer = 	new Observer(system, "Observer", 
										observer_particles, 0, ParticlesGeneratorOps::direct);
	/**
	 * @brief 	Body contact map.
	 * @details The contact map gives the data conntections between the bodies.
	 * 			Basically the the rang of bidies to build neighbor particle lists.
	 */
	SPHBodyTopology body_topology = { { water_block, { wall_boundary, inserted_body } },
									  { wall_boundary, { } }, { inserted_body, { water_block } }, 
									  { beam_observer,{inserted_body} } };
	system.SetBodyTopology(&body_topology);
	/**
	 * @brief 	Simulation set up.
	 */
	system.SetupSPHSimulation();
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	/**
	 * @brief 	Methods used only once.
	 */
	/** Obtain the initial number density of fluid. */
	fluid_dynamics::InitialNumberDensity 		fluid_initial_number_density(water_block, { wall_boundary, inserted_body });
	/** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_boundary, {});
	/** Initialize normal direction of the inserted body. */
	solid_dynamics::NormalDirectionSummation 	get_inserted_body_normal(inserted_body, {});
	/** Corrected strong configuration. */
	solid_dynamics::CorrectConfiguration 		inserted_body_corrected_configuration_in_strong_form(inserted_body, {});
	/**
	 * @brief 	Methods used for time stepping.
	 */
	/** Initialize particle acceleration. */
	InitializeOtherAccelerations 	initialize_other_acceleration(water_block);
	/** Periodic bounding. */
	PeriodicBoundingInXDirection 	periodic_bounding(water_block);
	/** Periodic BCs. */
	PeriodicConditionInXDirection 	periodic_condition(water_block);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_desnity(water_block, { wall_boundary, inserted_body });
	/** Divergence correction. */
	fluid_dynamics::DivergenceCorrection 		divergence_correction(water_block, { wall_boundary });
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize 	get_fluid_adevction_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationVerlet 		pressure_relaxation(water_block, { wall_boundary, inserted_body });
	/** Computing viscous acceleration. */
	fluid_dynamics::ComputingViscousAcceleration 	viscous_acceleration(water_block, { wall_boundary, inserted_body });
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrection 	transport_velocity_correction(water_block, { wall_boundary, inserted_body });
	/** Computing vorticity in the flow. */
	fluid_dynamics::ComputingVorticityInFluidField 	compute_vorticity(water_block);
		/** Inflow boundary condition. */
	ParabolicInflow 								
		parabolic_inflow(water_block, new InflowBuffer(water_block, "Buffer"));
	/** 
	 * @brief Algorithms of FSI. 
	 */
	/** Compute the force exerted on elastic gate due to fluid pressure and visocisty. */
	solid_dynamics::FluidPressureForceOnSolid 	fluid_pressure_force_on_insrted_body(inserted_body, { water_block }, &fluid);
	solid_dynamics::FluidViscousForceOnSolid 	fluid_viscous_force_on_insrted_body(inserted_body, { water_block }, &fluid);
	/** 
	 * @brief Algorithms of Elastic dynamics. 
	 */
	/** Compute time step size of elastic solid. */
	solid_dynamics::GetAcousticTimeStepSize 	inserted_body_computing_time_step_size(inserted_body);
	/** Stress relaxation for the beam. */
	solid_dynamics::StressRelaxationFirstStep 	inserted_body_stress_relaxation_first_step(inserted_body);
	solid_dynamics::StressRelaxationSecondStep 	inserted_body_stress_relaxation_second_step(inserted_body);
	/** Constrain region of the inserted body. */
	solid_dynamics::ConstrainSolidBodyRegion
		constrain_beam_base(inserted_body, new BeamBase(inserted_body, "BeamBase"));
	/** Computing the verage velocity. */
	solid_dynamics::InitializeDisplacement 			inserted_body_initialize_displacement(inserted_body);
	solid_dynamics::UpdateAverageVelocity 			inserted_body_average_velocity(inserted_body);
	/** Update norm .*/
	solid_dynamics::UpdateElasticNormalDirection 	inserted_body_update_normal(inserted_body);
	/**
	 * @brief 	Methods used for updating data structure.
	 */
	/** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList 			update_water_block_cell_linked_list(water_block);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsConfiguration 			update_water_block_configuration(water_block);
	/** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList 			update_inserted_body_cell_linked_list(inserted_body);
	/** Update the contact configuration for a given contact map. */
	ParticleDynamicsContactConfiguration 	update_inserted_body_contact_configuration(inserted_body);
	/**
	 * @brief The time stepping starts here.
	 */
	/** Pre-simultion*/
	periodic_condition.parallel_exec();
	update_water_block_configuration.parallel_exec();
	get_wall_normal.parallel_exec();
	get_inserted_body_normal.parallel_exec();
	fluid_initial_number_density.parallel_exec();
	inserted_body_corrected_configuration_in_strong_form.parallel_exec();
	/**
	 * @brief Output.
	 */
	Output 								output(system);
	WriteBodyStatesToVtu 				write_real_body_states_to_vtu(output, system.real_bodies_);
	WriteBodyStatesToPlt 				write_real_body_states_to_plt(output, system.real_bodies_);
	WriteTotalViscousForceOnSolid 		write_total_viscous_force_on_inserted_body(output, inserted_body);
	WriteObservedElasticDisplacement 	write_beam_tip_displacement(output, beam_observer, { inserted_body });
	WriteRestartFileToXml				write_restart_body_states(output, system.real_bodies_);
	/**
	 * @brief The time stepping starts here.
	 */
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if(system.restart_step_ != 0)
	{
		system.ResetSPHSimulationFromRestart();	

		update_water_block_cell_linked_list.parallel_exec();
		update_inserted_body_cell_linked_list.parallel_exec();

		periodic_condition.parallel_exec();
		update_water_block_configuration.parallel_exec();
		update_inserted_body_contact_configuration.parallel_exec();

		get_inserted_body_normal.parallel_exec();	
	}
	write_real_body_states_to_plt.WriteToFile(GlobalStaticVariables::physical_time_);
	write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	int rst_out = 5000;				/**< time step frquency for restart output. */
	Real End_Time = 200.0;			/**< End time. */
	Real D_Time = End_Time/200.0;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes. */
	Real dt = 0.0; 					/**< Default accoustic time step sizes. */
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Output global basic parameters. */
	output.WriteCaseSetup(GlobalStaticVariables::physical_time_, D_Time, End_Time);
	/**
	 * @brief Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integeral_time < D_Time) {

			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			divergence_correction.parallel_exec();
			initialize_other_acceleration.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			/** FSI for viscous force. */
			fluid_viscous_force_on_insrted_body.parallel_exec();
			/** Update normal direction on elastic body.*/
			inserted_body_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
					if (ite % rst_out == 0)
						write_restart_body_states.WriteToFile(Real(ite));
				}
				/** Fluid pressure relaxation. */
				pressure_relaxation.parallel_exec(dt);
				/** FSI for pressure force. */
				fluid_pressure_force_on_insrted_body.parallel_exec();
				/** Solid dynamics. */
				Real dt_s_sum = 0.0;
				inserted_body_initialize_displacement.parallel_exec();
				while (dt_s_sum < dt) {

					Real dt_s = inserted_body_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

					if (ite % 100 == 0) {
						cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_ << "	dt_s: "
							<< dt_s << "\n";
					}
					inserted_body_stress_relaxation_first_step.parallel_exec(dt_s);
					constrain_beam_base.parallel_exec();
					inserted_body_stress_relaxation_second_step.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				inserted_body_average_velocity.parallel_exec(dt);

				ite++;
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();
			}
			/** Water block confifuration and periodic constion. */
			periodic_bounding.parallel_exec();
			update_water_block_cell_linked_list.parallel_exec();
			periodic_condition.parallel_exec();
			update_water_block_configuration.parallel_exec();
			/** Inserted body contact configuration. */
			update_inserted_body_cell_linked_list.parallel_exec();
			update_inserted_body_contact_configuration.parallel_exec();
			write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
		}
		compute_vorticity.parallel_exec();
		tick_count t2 = tick_count::now();
		write_real_body_states_to_vtu.WriteToFile(GlobalStaticVariables::physical_time_);
		write_total_viscous_force_on_inserted_body.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
