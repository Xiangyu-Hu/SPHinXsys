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
//----------------------------------------------------------------------
//	Global parameters on the material properties
//----------------------------------------------------------------------
Real diffusion_coff = 1.0e-3;
Real bias_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
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
	//geometry
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
/** create a inflow buffer shape. */
MultiPolygon createInflowBufferShape()
{
	std::vector<Vecd> inflow_buffer_shape;
	inflow_buffer_shape.push_back(Vecd(0.0 - DL_sponge, 0.0));
	inflow_buffer_shape.push_back(Vecd(0.0 - DL_sponge, DH));
	inflow_buffer_shape.push_back(Vecd(0, DH));
	inflow_buffer_shape.push_back(Vecd(0, 0));
	inflow_buffer_shape.push_back(Vecd(0.0 - DL_sponge, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
/**  Thermo fluid body definition */
class ThermofluidBody : public FluidBody
{
public:
	ThermofluidBody(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/**  Thermo solid body definition */
class ThermosolidBody : public SolidBody
{
public:
	ThermosolidBody(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion fluid body
//----------------------------------------------------------------------
class ThermofluidBodyMaterial
	: public DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>
{
public:
	ThermofluidBodyMaterial()
		: DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>({"Phi"}, rho0_f, c_f, mu_f)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
	};
};
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion solid body
//----------------------------------------------------------------------
class ThermosolidBodyMaterial
	: public DiffusionReaction<SolidParticles, Solid>
{
public:
	ThermosolidBodyMaterial()
		: DiffusionReaction<SolidParticles, Solid>({"Phi"})
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
	};
};
//----------------------------------------------------------------------
//	Application dependent solid body initial condition
//----------------------------------------------------------------------
class ThermosolidBodyInitialCondition
	: public DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{

		if (-BW <= pos_n_[index_i][1] && pos_n_[index_i][1] <= 0.0)
		{
			species_n_[phi_][index_i] = phi_lower_wall;
		}

		if (DH <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH + BW)
		{
			species_n_[phi_][index_i] = phi_upper_wall;
		}
	};

public:
	explicit ThermosolidBodyInitialCondition(SolidBody &diffusion_solid_body)
		: DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>(diffusion_solid_body)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};
//----------------------------------------------------------------------
//	Application dependent fluid body initial condition
//----------------------------------------------------------------------
class ThermofluidBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidBody, FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{

		if (0 <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH)
		{
			species_n_[phi_][index_i] = phi_fluid_initial;
		}
	};

public:
	explicit ThermofluidBodyInitialCondition(FluidBody &diffusion_fluid_body)
		: DiffusionReactionInitialCondition<FluidBody, FluidParticles, WeaklyCompressibleFluid>(diffusion_fluid_body)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
class ThermalRelaxationComplex
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  FluidBody, FluidParticles, WeaklyCompressibleFluid,
		  RelaxationOfAllDiffussionSpeciesComplex<
			  FluidBody, FluidParticles, WeaklyCompressibleFluid, SolidBody, SolidParticles, Solid>,
		  ComplexBodyRelation>
{
public:
	explicit ThermalRelaxationComplex(ComplexBodyRelation &body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation){};
	virtual ~ThermalRelaxationComplex(){};
};
//----------------------------------------------------------------------
//	Case dependent inflow boundary condition.
//----------------------------------------------------------------------
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(FluidBody &fluid_body, BodyPartByCell &constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region),
		  u_ave_(0.0), u_ref_(1.0), t_ref(2.0) {}

	Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0)
		{
			u = 6.0 * u_ave_ * position[1] * (DH - position[1]) / DH / DH;
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
//	An observer particle generator
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** A measuring point at the center of the channel */
		Vec2d point_coordinate(0.0, DH * 0.5);
		positions_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
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
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	ThermofluidBody thermofluid_body(system, "ThermofluidBody");
	DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>
		diffusion_fluid_body_particles(thermofluid_body, makeShared<ThermofluidBodyMaterial>());

	ThermosolidBody thermosolid_body(system, "ThermosolidBody");
	DiffusionReactionParticles<SolidParticles, Solid>
		diffusion_solid_body_particles(thermosolid_body, makeShared<ThermosolidBodyMaterial>());

	ObserverBody temperature_observer(system, "FluidObserver");
	ObserverParticles temperature_observer_particles(temperature_observer, makeShared<ObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner fluid_body_inner(thermofluid_body);
	BodyRelationInner solid_body_inner(thermosolid_body);
	ComplexBodyRelation fluid_body_complex(fluid_body_inner, {&thermosolid_body});
	BodyRelationContact fluid_observer_contact(temperature_observer, {&thermofluid_body});

	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition(thermofluid_body, xAxis);
	ThermosolidBodyInitialCondition thermosolid_condition(thermosolid_body);
	ThermofluidBodyInitialCondition thermofluid_initial_condition(thermofluid_body);
	/** Corrected configuration for diffusion solid body. */
	solid_dynamics::CorrectConfiguration correct_configuration(solid_body_inner);
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(thermofluid_body);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationComplex update_density_by_summation(fluid_body_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step(thermofluid_body, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step(thermofluid_body);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<FluidBody, FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step(thermofluid_body);
	/** Diffusion process between two diffusion bodies. */
	ThermalRelaxationComplex thermal_relaxation_complex(fluid_body_complex);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(fluid_body_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(fluid_body_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(fluid_body_complex);
	/** Apply transport velocity formulation. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(fluid_body_complex);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInner compute_vorticity(fluid_body_inner);
	/** Inflow boundary condition. */
	MultiPolygonShape inflow_buffer_shape(createInflowBufferShape());
	BodyRegionByCell inflow_buffer(thermofluid_body, "Buffer", inflow_buffer_shape);
	ParabolicInflow parabolic_inflow(thermofluid_body, inflow_buffer);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states(in_output, system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<Real>>
		write_fluid_phi("Phi", in_output, fluid_observer_contact);
	ObservedQuantityRecording<Vecd>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);
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
	diffusion_solid_body_particles.initializeNormalDirectionFromBodyShape();
	correct_configuration.parallel_exec();
	thermosolid_condition.parallel_exec();
	thermofluid_initial_condition.parallel_exec();
	Real dt_thermal = get_thermal_time_step.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	Real End_Time = 10;
	Real D_Time = End_Time / 100.0; /**< time stamps for output,WriteToFile*/
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
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

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
			thermofluid_body.updateCellLinkedList();
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
