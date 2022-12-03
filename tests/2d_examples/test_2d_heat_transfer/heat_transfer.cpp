/**
 * @file 	heat_transfer.cpp
 * @brief 	This is a test to validate heat transfer between a flow and channel walls.
 * @author 	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;							/**< Channel length. */
Real DH = 0.4;							/**< Channel height. */
Real resolution_ref = DH / 25.0;		/**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
/** Boundary width, determined by specific layer of boundary particles. */
Real BW = resolution_ref * 4.0; /** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
// temperature observer location
StdVec<Vecd> observation_location = {Vecd(0.0, DH * 0.5)};
//----------------------------------------------------------------------
//	Global parameters on the material properties
//----------------------------------------------------------------------
Real diffusion_coff = 1.0e-3;
Real rho0_f = 1.0;					/**< Density. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Speed of sound. */
Real Re = 100.0;					/**< Reynolds number100. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the initial condition
//----------------------------------------------------------------------
Real phi_upper_wall = 20.0;
Real phi_lower_wall = 40.0;
Real phi_fluid_initial = 20.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createShape()
{
	// geometry
	std::vector<Vecd> shape;
	shape.push_back(Vecd(0.0 - DL_sponge, 0.0));
	shape.push_back(Vecd(0.0 - DL_sponge, DH));
	shape.push_back(Vecd(DL, DH));
	shape.push_back(Vecd(DL, 0.0));
	shape.push_back(Vecd(0.0 - DL_sponge, 0.0));
	return shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
	return outer_wall_shape;
}
/** create inner wall shape */
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

	return inner_wall_shape;
}
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
//----------------------------------------------------------------------
//	Case-dependent geometries
//----------------------------------------------------------------------
class ThermofluidBody : public MultiPolygonShape
{
public:
	explicit ThermofluidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createShape(), ShapeBooleanOps::add);
	}
};
class ThermosolidBody : public MultiPolygonShape
{
public:
	explicit ThermosolidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion fluid body
//----------------------------------------------------------------------
class ThermofluidBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
public:
	ThermofluidBodyMaterial()
		: DiffusionReaction<WeaklyCompressibleFluid>({"Phi"}, rho0_f, c_f, mu_f)
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
	};
};
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion solid body
//----------------------------------------------------------------------
class ThermosolidBodyMaterial : public DiffusionReaction<Solid>
{
public:
	ThermosolidBodyMaterial() : DiffusionReaction<Solid>({"Phi"})
	{
		// only default property is given, as no heat transfer within solid considered here.
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi");
	};
};
//----------------------------------------------------------------------
//	Application dependent solid body initial condition
//----------------------------------------------------------------------
class ThermosolidBodyInitialCondition
	: public DiffusionReactionInitialCondition<SolidParticles, Solid>
{
protected:
	size_t phi_;

public:
	explicit ThermosolidBodyInitialCondition(SPHBody &sph_body)
		: DiffusionReactionInitialCondition<SolidParticles, Solid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		if (-BW <= pos_[index_i][1] && pos_[index_i][1] <= 0.0)
		{
			species_n_[phi_][index_i] = phi_lower_wall;
		}

		if (DH <= pos_[index_i][1] && pos_[index_i][1] <= DH + BW)
		{
			species_n_[phi_][index_i] = phi_upper_wall;
		}
	};
};
//----------------------------------------------------------------------
//	Application dependent fluid body initial condition
//----------------------------------------------------------------------
class ThermofluidBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

public:
	explicit ThermofluidBodyInitialCondition(SPHBody &sph_body)
		: DiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		if (0 <= pos_[index_i][1] && pos_[index_i][1] <= DH)
		{
			species_n_[phi_][index_i] = phi_fluid_initial;
		}
	};
};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
class ThermalRelaxationComplex
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  RelaxationOfAllDiffusionSpeciesComplex<FluidParticles, WeaklyCompressibleFluid, SolidParticles, Solid>>
{
public:
	explicit ThermalRelaxationComplex(ComplexRelation &body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation){};
	virtual ~ThermalRelaxationComplex(){};
};
//----------------------------------------------------------------------
//	Case dependent inflow boundary condition.
//----------------------------------------------------------------------
class ParabolicInflow : public fluid_dynamics::InflowVelocityCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(BodyAlignedBoxByCell &aligned_box_part)
		: InflowVelocityCondition(aligned_box_part),
		  u_ave_(0.0), u_ref_(1.0), t_ref(2.0) {}

	Vecd getPrescribedVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0)
		{
			u = 1.5 * u_ave_ * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
			v = 0.0;
		}
		return Vecd(u, v);
	}
	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody thermofluid_body(system, makeShared<ThermofluidBody>("ThermofluidBody"));
	thermofluid_body.defineParticlesAndMaterial<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermofluidBodyMaterial>();
	thermofluid_body.generateParticles<ParticleGeneratorLattice>();

	SolidBody thermosolid_body(system, makeShared<ThermosolidBody>("ThermosolidBody"));
	thermosolid_body.defineParticlesAndMaterial<DiffusionReactionParticles<SolidParticles, Solid>, ThermosolidBodyMaterial>();
	thermosolid_body.generateParticles<ParticleGeneratorLattice>();

	ObserverBody temperature_observer(system, "FluidObserver");
	temperature_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation fluid_body_inner(thermofluid_body);
	InnerRelation solid_body_inner(thermosolid_body);
	ComplexRelation fluid_body_complex(fluid_body_inner, {&thermosolid_body});
	ContactRelation fluid_observer_contact(temperature_observer, {&thermofluid_body});

	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	PeriodicConditionUsingCellLinkedList periodic_condition(thermofluid_body, thermofluid_body.getBodyShapeBounds(), xAxis);
	SimpleDynamics<ThermosolidBodyInitialCondition> thermosolid_condition(thermosolid_body);
	SimpleDynamics<ThermofluidBodyInitialCondition> thermofluid_initial_condition(thermofluid_body);
	SimpleDynamics<NormalDirectionFromBodyShape> thermosolid_body_normal_direction(thermosolid_body);
	/** Initialize particle acceleration. */
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(thermofluid_body);
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_body_complex);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step(thermofluid_body, U_f);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step(thermofluid_body);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step(thermofluid_body);
	/** Diffusion process between two diffusion bodies. */
	ThermalRelaxationComplex thermal_relaxation_complex(fluid_body_complex);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(fluid_body_complex);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(fluid_body_complex);
	/** Computing viscous acceleration. */
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(fluid_body_complex);
	/** Apply transport velocity formulation. */
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transport_velocity_correction(fluid_body_complex);
	/** Computing vorticity in the flow. */
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_body_inner);
	/** Inflow boundary condition. */
	BodyAlignedBoxByCell inflow_buffer(
		thermofluid_body, makeShared<AlignedBoxShape>(Transform2d(Vec2d(buffer_translation)), buffer_halfsize));
	SimpleDynamics<ParabolicInflow, BodyAlignedBoxByCell> parabolic_inflow(inflow_buffer);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<Real>>
		write_fluid_phi("Phi", io_environment, fluid_observer_contact);
	ObservedQuantityRecording<Vecd>
		write_fluid_velocity("Velocity", io_environment, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
	 * but before the configuration build up. */
	periodic_condition.update_cell_linked_list_.parallel_exec();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the wall. */
	thermosolid_body_normal_direction.parallel_exec();
	thermosolid_condition.parallel_exec();
	thermofluid_initial_condition.parallel_exec();
	Real dt_thermal = get_thermal_time_step.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	Real end_time = 10;
	Real output_interval = end_time / 100.0; /**< time stamps for output,WriteToFile*/
	int number_of_iterations = 0;
	int screen_output_interval = 40;
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
			Real Dt = get_fluid_advection_time_step.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec();

			size_t inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt = SMIN(SMIN(dt_thermal, get_fluid_time_step.parallel_exec()), Dt);
				pressure_relaxation.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				thermal_relaxation_complex.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_condition.bounding_.parallel_exec();
			thermofluid_body.updateCellLinkedListWithParticleSort(100);
			periodic_condition.update_cell_linked_list_.parallel_exec();
			fluid_body_complex.updateConfiguration();
		}
		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		fluid_observer_contact.updateConfiguration();
		write_real_body_states.writeToFile();
		write_fluid_phi.writeToFile(number_of_iterations);
		write_fluid_velocity.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}

	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_fluid_phi.newResultTest();

	return 0;
}
