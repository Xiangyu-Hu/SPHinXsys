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
Real BW = particle_spacing_ref * 4.0;       /**< Extending width for BCs. */

Real cx = 2.0;                              /**< Center of fish in x direction. */
Real cy = 4.0;                              /**< Center of fish in y direction. */
Real fish_length = 3.738;                   /**< Length of fish. */
Vec3d tethering_point(-1.0, cy, 0.0);        /**< The tethering point. */
											/**
											* Material properties of the fluid.
											*/
Real rho0_f = 1.0;
Real U_f = 1.0;
Real c_f = 10.0*U_f;
Real Re = 5.0e3;
Real mu_f = rho0_f * U_f * (fish_length) / Re;
Real k_f = 0.0;
Real initial_pressure = 0.0;
Vec2d intial_velocity(0.0, 0.0);
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
	pnts_shaping_inner_wall.push_back(Point(-DLsponge - 2.0*BW, 0.0));
	pnts_shaping_inner_wall.push_back(Point(-DLsponge - 2.0*BW, DH));
	pnts_shaping_inner_wall.push_back(Point(DL + 2.0*BW, DH));
	pnts_shaping_inner_wall.push_back(Point(DL + 2.0*BW, 0.0));
	pnts_shaping_inner_wall.push_back(Point(-DLsponge - 2.0*BW, 0.0));

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
		body_region_.add_geometry(new Geometry(CreatWaterBlockShape()), RegionBooleanOps::add);
		/** Exclude the fish body. */
		std::vector<Point> fish_shape = CreatFishShape(cx, cy, fish_length, particle_spacing_*0.5);
		body_region_.add_geometry(new Geometry(fish_shape), RegionBooleanOps::sub);
		body_region_.done_modeling();
	}

	void InitialCondition()
	{
		SetAllParticleAtRest();
	}
};
/**
* Solid wall.
*/
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name,
		SolidBodyParticles &solid_particles, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, solid_particles,
			refinement_level, op)
	{
		body_region_.add_geometry(new Geometry(CreatOuterWallShape()), RegionBooleanOps::add);
		body_region_.add_geometry(new Geometry(CreatInnerWallShape()), RegionBooleanOps::sub);
		body_region_.done_modeling();

	}

	void InitialCondition()
	{
		SetAllParticleAtRest();
	}
};
/**
* Fish body with tethering constraint.
*/
class FishBody : public ElasticBody
{

public:
	FishBody(SPHSystem &system, string body_name, ElasticSolid* material,
		ElasticBodyParticles &elastic_particles, int refinement_level, ParticlesGeneratorOps op)
		: ElasticBody(system, body_name, material, elastic_particles,
			refinement_level, op)
	{

		std::vector<Point> fish_shape = CreatFishShape(cx, cy, fish_length, particle_spacing_);
		body_region_.add_geometry(new Geometry(fish_shape), RegionBooleanOps::add);
		body_region_.done_modeling();
	}

	void InitialCondition()
	{
		SetAllParticleAtRest();
	}
};
/**
* @brief create fish head for constraint
*/
class FishHead : public SolidBodyPartForSimbody
{
public:
	FishHead(SolidBody *solid_body, 
		string constrianed_region_name, Real solid_body_density)
		: SolidBodyPartForSimbody(solid_body, 
			constrianed_region_name, solid_body_density)
	{
		//geometry
		std::vector<Point> fish_shape = CreatFishShape(cx, cy, fish_length, solid_body_->particle_spacing_);
		soild_body_part_region_.add_geometry(new Geometry(fish_shape), RegionBooleanOps::add);
		soild_body_part_region_.add_geometry(new Geometry(CreatFishBlockingShape()), RegionBooleanOps::sub);
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
* Definition of an observer body with serveral particle to record particular property.
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
	ParabolicInflow(WeaklyCompressibleFluidBody* fluid_body,
		FluidBodyPart *constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region)
	{
		u_ave_ = 0.0;
		u_ref_ = 1.0;
		t_ref = 4.0;
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
		u_ave_ = run_time < t_ref ? 0.5*u_ref_*(1.0 - cos(pi*run_time / t_ref)) : u_ref_;
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
	SPHSystem system(Vec2d(-DLsponge - BW, -BW),
		Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/**
	* @brief   Particles and body creation for water.
	*/
	SymmetricTaitFluid                  fluid("Water", rho0_f, c_f, mu_f, k_f);
	WeaklyCompressibleFluidParticles    fluid_particles("WaterBody");
	WaterBlock *water_block = new WaterBlock(system, "WaterBody",
		&fluid, fluid_particles, 0, ParticlesGeneratorOps::lattice);
	/**
	* @brief   Particles and body creation for wall boundary.
	*/
	SolidBodyParticles                  solid_particles("Wall");
	WallBoundary *wall_boundary = new   WallBoundary(system, "Wall",
		solid_particles, 0, ParticlesGeneratorOps::lattice);
	/**
	* @brief   Particles and body creation for fish.
	*/
	NeoHookeanSolid             solid_material("NeoHookeanSolid",
		rho0_s, Youngs_modulus, poisson, 0.0);
	ElasticBodyParticles        fish_body_particles("FishBody");
	FishBody *fish_body = new   FishBody(system, "FishBody", &solid_material,
		fish_body_particles, 1, ParticlesGeneratorOps::lattice);
	/**
	* @brief   Particle and body creation of gate observer.
	*/
	ObserverParticles           observer_particles("Observer");
	Observer *fish_observer = new Observer(system, "Observer",
		observer_particles, 0, ParticlesGeneratorOps::direct);
	/**
	* @brief   Body contact map.
	* @details The contact map gives the data conntections between the bodies.
	*          Basically the the range of bidies to build neighbor particle lists.
	*/
	SPHBodyTopology body_topology = { { water_block,{ wall_boundary, fish_body } },
	{ wall_boundary,{} },{ fish_body,{ water_block } },{ fish_observer,{ fish_body } } };
	system.SetBodyTopology(&body_topology);
	/**
	* Setting up the simulation.
	*/
	system.SetupSPHSimulation();
	/**
	* This section define all numerical methods will be used in this case.
	*/
	/** Methods only used only once.*/
	fluid_dynamics::InitialNumberDensity
		fluid_initial_number_density(water_block, { wall_boundary, fish_body });
	/** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation get_wall_normal(wall_boundary, {});
	/** Initialize normal direction of the inserted body. */
	solid_dynamics::NormalDirectionReNormalization get_fish_body_normal(fish_body, {});
	/** Corrected strong configuration.*/
	solid_dynamics::CorrectConfiguration
		fish_body_corrected_configuration_in_strong_form(fish_body, {});
	/**
	* Common paritcle dynamics.
	*/
	InitializeOtherAccelerations
		initialize_fluid_acceleration(water_block);
	/** Periodic bounding. */
	PeriodicBoundingInXDirection
		periodic_bounding(water_block);
	/** Periodic boundary condition. */
	PeriodicConditionInXDirection
		periodic_condition(water_block);

	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation
		update_fluid_desnity(water_block, { wall_boundary, fish_body });
	//divergence correction
	fluid_dynamics::DivergenceCorrection
		divergence_correction(water_block, { wall_boundary });
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationVerlet
		pressure_relaxation(water_block, { wall_boundary, fish_body });
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
		fluid_pressure_force_on_fish_body(fish_body, { water_block }, &fluid);
	solid_dynamics::FluidViscousForceOnSolid
		fluid_viscous_force_on_fish_body(fish_body, { water_block }, &fluid);
	/**
	* Solid dynmaics.
	*/
	/** Time step size caclutation. */
	solid_dynamics::GetAcousticTimeStepSize fish_body_computing_time_step_size(fish_body);
	/** Process of stress relaxation. */
	solid_dynamics::StressRelaxationFirstStep
		fish_body_stress_relaxation_first_step(fish_body);
	solid_dynamics::StressRelaxationSecondStep
		fish_body_stress_relaxation_second_step(fish_body);
	/** Update normal direction on fish body.*/
	solid_dynamics::UpdateElasticNormalDirection
		fish_body_update_normal(fish_body);
	/** Compute the average velocity on fish body. */
	solid_dynamics::InitializeDisplacement
		fish_body_initialize_displacement(fish_body);
	solid_dynamics::UpdateAverageVelocity
		fish_body_average_velocity(fish_body);
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
	FishHead *fish_head = new FishHead(fish_body, "FishHead", rho0_s);
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
	SimTK::Visualizer viz(MBsystem);
	MBsystem.addEventReporter(new Visualizer::Reporter(viz, 0.01));
	/** Initialize the system and state. */
	SimTK::State state = MBsystem.realizeTopology();
	viz.report(state);
	cout << "Hit ENTER to run a short simulation ...";
	getchar();
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

	/** Output. */
	Output output(system);
	WriteBodyStatesToVtu        write_real_body_states(output, system.real_bodies_);
	WriteTotalForceOnSolid      write_total_force_on_fish(output, fish_body);
	WriteObservedElasticDisplacement write_fish_displacement(output, fish_observer, { fish_body });
	/**
	* Time steeping starts here.
	*/
	GlobalStaticVariables::physical_time_ = 0.0;
	/**
	* Initial periodic boundary condition which copies the particle identifies
	* as extra cell linked list form periodic regions to the corresponding boundaries
	* for buiding up of extra configuration.
	*/
	periodic_condition.parallel_exec();
	/** Update configuration after periodic boundary condition.*/
	update_water_block_configuration.parallel_exec();
	/** Prepare quantities, e.g. wall normal, fish body norm,
	* fluid initial number density and configuration of fish particles, will be used once only.
	*/
	get_wall_normal.parallel_exec();
	get_fish_body_normal.parallel_exec();
	fluid_initial_number_density.parallel_exec();
	fish_body_corrected_configuration_in_strong_form.parallel_exec();
	/** Output for initial condition. */
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_fish_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	* Time parameters
	*/
	int ite = 0;
	Real End_Time = 200.0;
	Real D_Time = End_Time / 200.0;
	Real Dt = 0.0;      /**< Default advection time step sizes. */
	Real dt = 0.0;      /**< Default accoustic time step sizes. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Output global basic parameters. */
	output.WriteCaseSetup(End_Time, D_Time, GlobalStaticVariables::physical_time_);
	/**
	* Main loop starts here.
	*/
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < D_Time)
		{
			initialize_fluid_acceleration.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			divergence_correction.parallel_exec();
			initialize_fluid_acceleration.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);
			transport_velocity_stress.parallel_exec();
			/** Viscous force exerting on fish body. */
			fluid_viscous_force_on_fish_body.parallel_exec();
			/** Update normal direction on fish body. */
			fish_body_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "   dt: "
						<< dt << "\n";
				}
				/** Fluid dynamics process. */
				pressure_relaxation.parallel_exec(dt);
				/** Fluid pressure force exerting on fish. */
				fluid_pressure_force_on_fish_body.parallel_exec();
				/** Relax fish body by solid dynamics. */
				Real dt_s_sum = 0.0;
				fish_body_initialize_displacement.parallel_exec();
				while (dt_s_sum < dt)
				{
					Real dt_s = fish_body_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

					if (ite % 100 == 0) {
						cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_ << "   dt_s: "
							<< dt_s << "\n";
					}
					fish_body_stress_relaxation_first_step.parallel_exec(dt_s);
					SimTK::State &state_for_update = integ.updAdvancedState();
					force_on_bodies.clearAllBodyForces(state_for_update);
					force_on_bodies.setOneBodyForce(state_for_update, tethered_spot,
						force_on_tethered_spot.parallel_exec());
					integ.stepBy(dt_s);
					constriant_tethered_spot.parallel_exec();
					fish_body_stress_relaxation_second_step.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				fish_body_average_velocity.parallel_exec(dt);
				write_total_force_on_fish.WriteToFile(GlobalStaticVariables::physical_time_);

				ite++;
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();

			}

			const State& s = integ.getState();
			viz.report(s);
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