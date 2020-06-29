/**
* @file    fish_passive_flapping_in_flow.cpp
* @brief   fish flapping passively in flow
* @author  Xiangyu Hu and Chi Zhang
* @version 0.1
*/
#include "sphinxsys.h"
/**
* Create the shapes for fish and bones.
*/
#include "fish_and_bones.h"
/**
* @brief Namespace cite here.
*/
using namespace SPH;
/**
* @brief Basic geometry parameters and numerical setup.
*/
Real DL = 11.0;                             /**< Channel length. */
Real DH = 8.0;                              /**< Channel height. */

Real particle_spacing_ref = 0.05;           /** Intial particle spacing. */
Real DLsponge = particle_spacing_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0;        /**< Extending width for BCs. */

Real cx = 2.0;                              /**< Center of fish in x direction. */
Real cy = 4.0;                              /**< Center of fish in y direction. */
Real fish_length = 3.738;                   /**< Length of fish. */
Vec3d tethering_point(-1.0, cy, 0.0);        /**< The tethering point. */
/**
 * Material properties of the fluid.
 */
Real rho0_f = 1.0;
Real U_f = 1.0;
Real c_f = 10.0 * U_f;
Real Re = 5.0e3;
Real mu_f = rho0_f * U_f * (fish_length) / Re;
/**
 * Material properties of the fish body.
 */
Real rho0_s = 1.0;
Real poisson = 0.49;
Real Ae = 2.0e2;
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
/**
 * Basic geometries for construct SPH bodies.
 */
/**
* @brief create a water block shape
*/
std::vector<Point> CreatWaterBlockShape()
{
	std::vector<Point> pnts_shaping_water_block;
	pnts_shaping_water_block.push_back(Point(-DLsponge, 0.0));
	pnts_shaping_water_block.push_back(Point(-DLsponge, DH));
	pnts_shaping_water_block.push_back(Point(DL, DH));
	pnts_shaping_water_block.push_back(Point(DL, 0.0));
	pnts_shaping_water_block.push_back(Point(-DLsponge, 0.0));

	return pnts_shaping_water_block;
}
/**
* @brief create a buffer for water block shape
*/
std::vector<Point> CreatInflowBufferShape()
{
	std::vector<Point> pnts_buffer;
	pnts_buffer.push_back(Point(-DLsponge, 0.0));
	pnts_buffer.push_back(Point(-DLsponge, DH));
	pnts_buffer.push_back(Point(0.0, DH));
	pnts_buffer.push_back(Point(0.0, 0.0));
	pnts_buffer.push_back(Point(-DLsponge, 0.0));

	return pnts_buffer;
}
/**
* @brief create outer wall shape
*/
std::vector<Point> CreatOuterWallShape()
{
	std::vector<Point> pnts_shaping_outer_wall;
	pnts_shaping_outer_wall.push_back(Point(-DLsponge - BW, -BW));
	pnts_shaping_outer_wall.push_back(Point(-DLsponge - BW, DH + BW));
	pnts_shaping_outer_wall.push_back(Point(DL + BW, DH + BW));
	pnts_shaping_outer_wall.push_back(Point(DL + BW, -BW));
	pnts_shaping_outer_wall.push_back(Point(-DLsponge - BW, -BW));

	return pnts_shaping_outer_wall;
}
/**
* @brief create inner wall shape
*/
std::vector<Point> CreatInnerWallShape()
{
	std::vector<Point> pnts_shaping_inner_wall;
	pnts_shaping_inner_wall.push_back(Point(-DLsponge - 2.0 * BW, 0.0));
	pnts_shaping_inner_wall.push_back(Point(-DLsponge - 2.0 * BW, DH));
	pnts_shaping_inner_wall.push_back(Point(DL + 2.0 * BW, DH));
	pnts_shaping_inner_wall.push_back(Point(DL + 2.0 * BW, 0.0));
	pnts_shaping_inner_wall.push_back(Point(-DLsponge - 2.0 * BW, 0.0));

	return pnts_shaping_inner_wall;
}
/**
* @brief create blocking shape to separate fish head out
*/
Real head_size = 1.0;
std::vector<Point> CreatFishBlockingShape()
{
	std::vector<Point> pnts_blocking_shape;
	pnts_blocking_shape.push_back(Point(cx + head_size, cy - 0.4));
	pnts_blocking_shape.push_back(Point(cx + head_size, cy + 0.4));
	pnts_blocking_shape.push_back(Point(cx + 5.0, cy + 0.4));
	pnts_blocking_shape.push_back(Point(cx + 5.0, cy - 0.4));
	pnts_blocking_shape.push_back(Point(cx + head_size, cy - 0.4));

	return pnts_blocking_shape;
}
/**
* Water body defintion.
*/
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: FluidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> water_block_shape = CreatWaterBlockShape();
		body_region_.add_geometry(new Geometry(water_block_shape), RegionBooleanOps::add);
		/** Exclude the fish body. */
		std::vector<Point> fish_shape = CreatFishShape(cx, cy, fish_length, particle_spacing_ * 0.5);
		body_region_.add_geometry(new Geometry(fish_shape), RegionBooleanOps::sub);
		body_region_.done_modeling();
	}
};
/**
 * @brief 	Case dependent material properties definition.
 */
class WaterMaterial : public SymmetricTaitFluid
{
public:
	WaterMaterial() : SymmetricTaitFluid()
	{
		rho_0_ = rho0_f;
		c_0_ = c_f;
		mu_ = mu_f;

		assignDerivedMaterialParameters();
	}
};
/**
* Solid wall.
*/
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem& system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> outer_wall_shape = CreatOuterWallShape();
		std::vector<Point> inner_wall_shape = CreatInnerWallShape();
		body_region_.add_geometry(new Geometry(outer_wall_shape), RegionBooleanOps::add);
		body_region_.add_geometry(new Geometry(inner_wall_shape), RegionBooleanOps::sub);
		body_region_.done_modeling();
	}
};
/**
* Fish body with tethering constraint.
*/
class FishBody : public SolidBody
{

public:
	FishBody(SPHSystem& system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> fish_shape = CreatFishShape(cx, cy, fish_length, particle_spacing_);
		body_region_.add_geometry(new Geometry(fish_shape), RegionBooleanOps::add);
		body_region_.done_modeling();
	}
};
/**
*@brief Define gate material.
*/
class FishMaterial : public NeoHookeanSolid
{
public:
	FishMaterial() : NeoHookeanSolid()
	{
		rho_0_ = rho0_s;
		E_0_ = Youngs_modulus;
		nu_ = poisson;

		assignDerivedMaterialParameters();
	}
};

/**
* @brief create fish head for constraint
*/
class FishHead : public SolidBodyPartForSimbody
{
public:
	FishHead(SolidBody* solid_body,
		string constrianed_region_name, Real solid_body_density)
		: SolidBodyPartForSimbody(solid_body,
			constrianed_region_name)
	{
		//geometry
		std::vector<Point> fish_shape = CreatFishShape(cx, cy, fish_length, body_->particle_spacing_);
		std::vector<Point> fish_blocking_shape = CreatFishBlockingShape();
		body_part_region_.add_geometry(new Geometry(fish_shape), RegionBooleanOps::add);
		body_part_region_.add_geometry(new Geometry(fish_blocking_shape), RegionBooleanOps::sub);
		//finish the region modeling
		body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
	}
};
/**
* @brief inflow buffer
*/
class InflowBuffer : public BodyPartByCell
{
public:
	InflowBuffer(FluidBody* fluid_body, string constrianed_region_name)
		: BodyPartByCell(fluid_body, constrianed_region_name)
	{
		/** Geomerty definition. */
		std::vector<Point> inflow_shape = CreatInflowBufferShape();
		body_part_region_.add_geometry(new Geometry(inflow_shape), RegionBooleanOps::add);
		/** Finalize the geometry definition and correspoding opertation. */
		body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartCells();
	}
};
/**
* Definition of an observer body with serveral particle to record particular property.
*/
class Observer : public FictitiousBody
{
public:
	Observer(SPHSystem& system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{
		/** postion and volume. */
		body_input_points_volumes_.push_back(make_pair(Point(cx + particle_spacing_ref, cy), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(cx + fish_length - particle_spacing_ref, cy), 0.0));
	}
};
/**
* Inflow boundary condition.
*/
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(FluidBody* fluid_body,
		BodyPartByCell* constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region)
	{
		u_ave_ = 0.0;
		u_ref_ = 1.0;
		t_ref = 4.0;
	}

	Vecd GetInflowVelocity(Vecd& position, Vecd& velocity)
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0) {
			u = 6.0 * u_ave_ * position[1] * (DH - position[1]) / DH / DH;
			v = 0.0;
		}
		return Vecd(u, v);
	}

	void PrepareConstraint() override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(pi * run_time / t_ref)) : u_ref_;
	}
};
/**
* Main program starts here.
*/
int main()
{
	/**
	* Build up context -- a SPHSystem.
	*/
	SPHSystem system(Vec2d(-DLsponge - BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;
	/**
	* @brief   Particles and body creation for water.
	*/
	WaterBlock* water_block = new WaterBlock(system, "WaterBody",
		0, ParticlesGeneratorOps::lattice);
	WaterMaterial    *water_fluid = new WaterMaterial();
	FluidParticles    fluid_particles(water_block, water_fluid);
	/**
	* @brief   Particles and body creation for wall boundary.
	*/
	WallBoundary* wall_boundary = new   WallBoundary(system, "Wall",
		0, ParticlesGeneratorOps::lattice);
	SolidParticles                  solid_particles(wall_boundary);
	/**
	* @brief   Particles and body creation for fish.
	*/
	FishBody* fish_body = new   FishBody(system, "FishBody", 1, ParticlesGeneratorOps::lattice);
	FishMaterial   *fish_body_material = new FishMaterial();
	ElasticSolidParticles  fish_body_particles(fish_body, fish_body_material);
	/**
	* @brief   Particle and body creation of gate observer.
	*/
	Observer* fish_observer = new Observer(system, "Observer",
		1, ParticlesGeneratorOps::direct);
	BaseParticles           observer_particles(fish_observer);
	/**
	* @brief   Body contact map.
	* @details The contact map gives the data conntections between the bodies.
	*          Basically the the range of bidies to build neighbor particle lists.
	*/
	SPHBodyTopology body_topology = { { water_block,{ wall_boundary, fish_body } },
	{ wall_boundary,{} },{ fish_body,{ water_block } },{ fish_observer,{ fish_body } } };
	system.SetBodyTopology(&body_topology);
	/** Output. */
	In_Output in_output(system);
	WriteBodyStatesToPlt        write_real_body_states(in_output, system.real_bodies_);
	WriteTotalForceOnSolid      write_total_force_on_fish(in_output, fish_body);
	WriteAnObservedQuantity<Vecd, BaseParticles,
		BaseParticleData, &BaseParticles::base_particle_data_, &BaseParticleData::pos_n_>
		write_fish_displacement("Displacement", in_output, fish_observer, fish_body);
	/**
	* @brief   Methods used for updating data structure.
	*/
	/** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList          update_water_block_cell_linked_list(water_block);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsConfiguration           update_water_block_configuration(water_block);
	/** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList          update_fish_body_cell_linked_list(fish_body);
	/** Update the contact configuration for a given contact map. */
	ParticleDynamicsContactConfiguration    update_fish_body_contact_configuration(fish_body);
	/** Update inner configuration of a body. */
	ParticleDynamicsInnerConfiguration update_fish_body_inner_configuration(fish_body);
	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(water_block, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(water_block, 0);

	/** check whether run particle relaxation for body fiited particle distribution. */
	if (system.run_particle_relaxation_) {
		/** add background level set for particle realxation. */
		fish_body->addBackgroundMesh();
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		 /** Random reset the insert body particle position. */
		RandomizePartilePosition  random_fish_body_particles(fish_body);
		/** Write backgroung level set. */
		WriteBodyMeshToPlt write_fish_body_background_mesh(in_output, fish_body);
		/** Write the body state to Vtu file. */
		WriteBodyStatesToPlt 		write_fish_body(in_output, { fish_body });
		/** Write the particle reload files. */
		WriteReloadParticle 		write_particle_reload_files(in_output, { fish_body });

		/** bounding particles to insert body surface. */
		relax_dynamics::BodySurfaceBounding
			body_surface_bounding(fish_body, new NearBodySurface(fish_body));
		/** Compute the time step for particle relaxation. */
		relax_dynamics::GetTimeStepSize get_solid_relax_timestep(fish_body);
		/** Physics relaxation algorithm. */
		relax_dynamics::PhysicsRelaxationInner 	relax_process_for_solid(fish_body);
		/** finilaizing  particle number desnity and inital position after relaxatoin. */
		relax_dynamics::FinalizingParticleRelaxation finalizing_inserted_body_particles(fish_body);
		/**
		  * @brief 	Particle relaxation starts here.
		  */
		write_fish_body_background_mesh.WriteToFile(0.0);
		update_fish_body_cell_linked_list.parallel_exec();
		update_fish_body_inner_configuration.parallel_exec();
		body_surface_bounding.parallel_exec();
		random_fish_body_particles.parallel_exec(0.25);
		write_fish_body.WriteToFile(0.0);

		/** relax particles of the insert body. */
		int ite_p = 0;
		Real dt_p = 0.0;
		while (ite_p < 1000)
		{
			dt_p = get_solid_relax_timestep.parallel_exec();

			relax_process_for_solid.parallel_exec(dt_p);
			body_surface_bounding.parallel_exec();
			ite_p += 1;

			update_fish_body_cell_linked_list.parallel_exec();
			update_fish_body_inner_configuration.parallel_exec();
			if (ite_p % 200 == 0)
			{
				cout << fixed << setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_fish_body.WriteToFile(Real(ite_p) * 1.0e-4);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
		finalizing_inserted_body_particles.parallel_exec();

		/** Output results. */
		write_particle_reload_files.WriteToFile(0.0);
		return 0;
	}

	/**
	* This section define all numerical methods will be used in this case.
	*/
	/** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation get_wall_normal(wall_boundary, {});
	/** Initialize normal direction of the inserted body. */
	solid_dynamics::NormalDirectionReNormalization get_fish_body_normal(fish_body, {});
	/** Corrected strong configuration.*/
	solid_dynamics::CorrectConfiguration
		fish_body_corrected_configuration_in_strong_form(fish_body);
	/**
	* Common paritcle dynamics.
	*/
	InitializeATimeStep
		initialize_a_fluid_step(water_block);

	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation
		update_fluid_desnity(water_block, { wall_boundary, fish_body });
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationFirstHalf
		pressure_relaxation_first_half(water_block, { wall_boundary, fish_body });
	fluid_dynamics::PressureRelaxationSecondHalf
		pressure_relaxation_second_half(water_block, { wall_boundary, fish_body });
	/** Computing viscous acceleration. */
	fluid_dynamics::ComputingViscousAcceleration
		viscous_acceleration(water_block, { wall_boundary, fish_body });
	/** Impose transport velocity formulation. */
	fluid_dynamics::TransportVelocityCorrection
		transport_velocity_correction(water_block, { wall_boundary, fish_body });
	fluid_dynamics::TransportVelocityStress
		transport_velocity_stress(water_block, { wall_boundary, fish_body });
	/** Computing vorticity in the flow. */
	fluid_dynamics::ComputingVorticityInFluidField
		compute_vorticity(water_block);
	/** Inflow boundary condition. */
	ParabolicInflow parabolic_inflow(water_block, new InflowBuffer(water_block, "Buffer"));

	/**
	* Fluid structure interaction model.
	*/
	solid_dynamics::FluidPressureForceOnSolid
		fluid_pressure_force_on_fish_body(fish_body, { water_block });
	solid_dynamics::FluidViscousForceOnSolid
		fluid_viscous_force_on_fish_body(fish_body, { water_block });
	/**
	* Solid dynmaics.
	*/
	/** Time step size caclutation. */
	solid_dynamics::GetAcousticTimeStepSize fish_body_computing_time_step_size(fish_body);
	/** Process of stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf
		fish_body_stress_relaxation_first_half(fish_body);
	solid_dynamics::StressRelaxationSecondHalf
		fish_body_stress_relaxation_second_half(fish_body);
	/** Update normal direction on fish body.*/
	solid_dynamics::UpdateElasticNormalDirection
		fish_body_update_normal(fish_body);
	/** Compute the average velocity on fish body. */
	solid_dynamics::InitializeDisplacement
		fish_body_initialize_displacement(fish_body);
	solid_dynamics::UpdateAverageVelocity
		fish_body_average_velocity(fish_body);
	/**
	* The multi body system from simbody.
	*/
	SimTK::MultibodySystem          MBsystem;
	/** The bodies or matter of the MBsystem. */
	SimTK::SimbodyMatterSubsystem   matter(MBsystem);
	/** The forces of the MBsystem.*/
	SimTK::GeneralForceSubsystem    forces(MBsystem);
	SimTK::CableTrackerSubsystem    cables(MBsystem);
	/** Mass proeprties of the fixed spot. */
	SimTK::Body::Rigid      fixed_spot_info(MassProperties(1.0, Vec3(0), UnitInertia(1)));
	FishHead* fish_head = new FishHead(fish_body, "FishHead", rho0_s);
	/** Mass properties of the consrained spot. */
	SimTK::Body::Rigid      tethered_spot_info(*fish_head->body_part_mass_properties_);
	/** Mobility of the fixed spot. */
	SimTK::MobilizedBody::Weld      fixed_spot(matter.Ground(), Transform(tethering_point),
		fixed_spot_info, Transform(Vec3(0)));
	/** Mobility of the tethered spot.
	  * Set the mass center as the origin location of the planar mobilizer
	  */
	Vec3 displacement0 = fish_head->initial_mass_center_ - tethering_point;
	SimTK::MobilizedBody::Planar    tethered_spot(fixed_spot,
		Transform(displacement0), tethered_spot_info, Transform(Vec3(0)));
	/** The tethering line give cable force.
	  * the start point of the cable path is at the origin loaction of the first mobilizer body,
	  * the end point is the tip of the fish head which has a distance to the origin
	  * location of the second mobilizer body origin location, here, the mass center
	  * of the fish head.
	  */
	Vec3 displacement_cable_end = Vec3(cx, cy, 0.0) - fish_head->initial_mass_center_;
	SimTK::CablePath    tethering_line(cables, fixed_spot,
		Vec3(0), tethered_spot, displacement_cable_end);
	SimTK::CableSpring  tethering_spring(forces, tethering_line, 100.0, 3.0, 10.0);

	//discreted forces acting on the bodies
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	fixed_spot_info.addDecoration(Transform(), DecorativeSphere(0.02));
	tethered_spot_info.addDecoration(Transform(), DecorativeSphere(0.4));
	/** Visulizer from simbody. */
	//SimTK::Visualizer viz(MBsystem);
	//MBsystem.addEventReporter(new Visualizer::Reporter(viz, 0.01));
	/** Initialize the system and state. */
	SimTK::State state = MBsystem.realizeTopology();
	//viz.report(state);
	/** Time steping method for multibody system.*/
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	integ.initialize(state);
	/**
	* Coupling between SimBody and SPH.
	*/
	solid_dynamics::ForceOnElasticBodyPartForSimBody
		force_on_tethered_spot(fish_body, fish_head,
			MBsystem, tethered_spot, force_on_bodies, integ);
	solid_dynamics::ConstrianSoildBodyPartBySimBody
		constriant_tethered_spot(fish_body,
			fish_head, MBsystem, tethered_spot, force_on_bodies, integ);

	/**
	* Time steeping starts here.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Using relaxed particle distribution if needed. */
	if (system.reload_particles_) {
		ReadReloadParticle		reload_insert_body_particles(in_output, { fish_body }, { "FishBody" });
		reload_insert_body_particles.ReadFromFile();
	}
	/**
	* Initial periodic boundary condition which copies the particle identifies
	* as extra cell linked list form periodic regions to the corresponding boundaries
	* for buiding up of extra configuration.
	*/
	system.InitializeSystemCellLinkedLists();
	periodic_condition.parallel_exec();
	system.InitializeSystemConfigurations();
	/** Prepare quantities, e.g. wall normal, fish body norm,
	* fluid initial number density and configuration of fish particles, will be used once only.
	*/
	get_wall_normal.parallel_exec();
	get_fish_body_normal.parallel_exec();
	fish_body_corrected_configuration_in_strong_form.parallel_exec();
	/** Output for initial condition. */
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_fish_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	* Time parameters
	*/
	int number_of_iterations = 0;
	int screen_output_interval = 100;
	Real End_Time = 200.0;
	Real D_Time = End_Time / 200.0;
	Real Dt = 0.0;      /**< Default advection time step sizes. */
	Real dt = 0.0;      /**< Default accoustic time step sizes. */
	Real dt_s = 0.0;	/**< Default acoustic time step sizes for solid. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	/**
	* Main loop starts here.
	*/
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);
			transport_velocity_stress.parallel_exec();
			/** Viscous force exerting on fish body. */
			fluid_viscous_force_on_fish_body.parallel_exec();
			/** Update normal direction on fish body. */
			fish_body_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				/** Fluid dynamics process, first half. */
				pressure_relaxation_first_half.parallel_exec(dt);
				/** Fluid pressure force exerting on fish. */
				fluid_pressure_force_on_fish_body.parallel_exec();
				/** Fluid dynamics process, second half. */
				pressure_relaxation_second_half.parallel_exec(dt);
				/** Relax fish body by solid dynamics. */
				Real dt_s_sum = 0.0;
				fish_body_initialize_displacement.parallel_exec();
				while (dt_s_sum < dt)
				{
					dt_s = fish_body_computing_time_step_size.parallel_exec();
					fish_body_stress_relaxation_first_half.parallel_exec(dt_s);
					SimTK::State& state_for_update = integ.updAdvancedState();
					force_on_bodies.clearAllBodyForces(state_for_update);
					force_on_bodies.setOneBodyForce(state_for_update, tethered_spot,
						force_on_tethered_spot.parallel_exec());
					integ.stepBy(dt_s);
					constriant_tethered_spot.parallel_exec();
					fish_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				fish_body_average_velocity.parallel_exec(dt);
				write_total_force_on_fish.WriteToFile(GlobalStaticVariables::physical_time_);

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();

			}
			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
			}
			number_of_iterations++;

			//const State& s = integ.getState();
			//viz.report(s);
			/** Water block confifuration and periodic constion. */
			periodic_bounding.parallel_exec();
			update_water_block_cell_linked_list.parallel_exec();
			periodic_condition.parallel_exec();
			update_water_block_configuration.parallel_exec();
			/** Fish body contact configuration. */
			update_fish_body_cell_linked_list.parallel_exec();
			update_fish_body_contact_configuration.parallel_exec();
			write_fish_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
		}
		tick_count t2 = tick_count::now();
		compute_vorticity.parallel_exec();
		write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_ * 0.001);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}