/**
 * @file 	heat_transfer.cpp
 * @brief 	This is a test to validate heat transfer between a flow and channel walls.
 * @author 	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 */
/** SPHinXsys Library. */
#include "sphinxsys.h"
using namespace SPH;


/** Geometry parameter. */
Real DL = 2.0;                /**< Channel length. */
Real DH = 0.4;              /**< Channel height. */
Real particle_spacing_ref = DH / 25.0;            /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 5.0;	/**< Sponge region to impose inflow condition. */

/**< Boundary width, determined by specific layer of boundary particles. */
Real BW = particle_spacing_ref * 4.0; 		

/** Material properties. */
Real diffusion_coff = 1.0e-3;
Real bias_diffusion_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));

/**
*@brief Material properties of the fluid.
*/
Real rho0_f = 1.0;		/**< Density. */
Real U_f =1.0;		 /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;	/**< Speed of sound. */
Real Re = 100.0;		/**< Reynolds number100. */
Real mu_f = rho0_f * U_f * DH / Re;	/**< Dynamics viscosity. */

/**
*@brief Temperatures.
*/
Real phi_upper_wall = 20.0;
Real phi_lower_wall = 40.0;
Real phi_fluid_initial = 20.0;

/** create a water block shape */
std::vector<Point> CreatShape()
{
	//geometry
	std::vector<Point> shape;
	shape.push_back(Point(0.0 - DL_sponge, 0.0));
	shape.push_back(Point(0.0 - DL_sponge, DH));
	shape.push_back(Point(DL, DH));
	shape.push_back(Point(DL, 0.0));
	shape.push_back(Point(0.0 - DL_sponge, 0.0));
	return shape;
}

/** create outer wall shape */
std::vector<Point> CreatOuterWallShape()
{
	std::vector<Point> outer_wall_shape;
	outer_wall_shape.push_back(Point(-DL_sponge - BW, -BW));
	outer_wall_shape.push_back(Point(-DL_sponge - BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, -BW));
	outer_wall_shape.push_back(Point(-DL_sponge - BW, -BW));
	return outer_wall_shape;
}

/** create inner wall shape */
std::vector<Point> CreatInnerWallShape()
{
	std::vector<Point> inner_wall_shape;
	inner_wall_shape.push_back(Point(-DL_sponge - 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Point(-DL_sponge - 2.0 * BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0 * BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Point(-DL_sponge - 2.0 * BW, 0.0));

	return inner_wall_shape;
}


/** create a water block buffer shape. */
std::vector<Point> CreatInflowBufferShape()
{
	std::vector<Point> inlfow_buffer_shape;
	inlfow_buffer_shape.push_back(Point(0.0 - DL_sponge, 0.0));
	inlfow_buffer_shape.push_back(Point(0.0 - DL_sponge, DH));
	inlfow_buffer_shape.push_back(Point(0, DH));
	inlfow_buffer_shape.push_back(Point(0, 0));
	inlfow_buffer_shape.push_back(Point(0.0 - DL_sponge, 0.0));

	return inlfow_buffer_shape;
}


/** inflow buffer */
class InflowBuffer : public BodyPartByCell
{
public:
	InflowBuffer(FluidBody* fluid_body, string constrained_region_name)
		: BodyPartByCell(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Point> inflow_buffer_shape = CreatInflowBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
		tagBodyPart();
	}
};



/**  Thermofluid body definition */
class ThermofluidBody : public FluidBody
{
public: 
	ThermofluidBody(SPHSystem &system, string body_name, int refinement_level)
		: FluidBody(system, body_name, refinement_level)
	{	
		std::vector<Point> body_shape = CreatShape();	
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(body_shape, ShapeBooleanOps::add);
	}
};

/**  Thermosolid body definition */
class ThermosolidBody : public SolidBody
{
public:
	ThermosolidBody(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		std::vector<Point>  outer_wall_shape = CreatOuterWallShape();
		std::vector<Point> inner_wall_shape = CreatInnerWallShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);

	}
};


/**
 * Setup heat condution material properties for diffusion fluid body 
 */
class ThermofluidBodyMaterial
	: public DiffusionReactionMaterial<FluidParticles, WeaklyCompressibleFluid>
{
public:
	ThermofluidBodyMaterial()
		: DiffusionReactionMaterial<FluidParticles, WeaklyCompressibleFluid>()
	{
		rho_0_ = rho0_f;
		c_0_ = c_f;
		mu_ = mu_f;

		//add a scalar for temperature in fluid
		insertASpecies("Phi");
		assignDerivedMaterialParameters();
		initializeDiffusion();
	}
	/** Initialize diffusion reaction material. */
	virtual void initializeDiffusion() override {
		DirectionalDiffusion* phi_diffusion
			= new DirectionalDiffusion(species_indexes_map_["Phi"], species_indexes_map_["Phi"],
				diffusion_coff, bias_diffusion_coff, bias_direction);
		species_diffusion_.push_back(phi_diffusion);
	};
};

/**
 * Setup heat conduction material properties for diffusion solid body
 */
class ThermosolidBodyMaterial
	: public DiffusionReactionMaterial<SolidParticles, Solid>
{
public:
	ThermosolidBodyMaterial()
		: DiffusionReactionMaterial<SolidParticles, Solid>()
	{
		//add a scalar for temperature in solid
		insertASpecies("Phi");
		assignDerivedMaterialParameters();
		initializeDiffusion();
	}
	/** Initialize diffusion reaction material. */
	virtual void initializeDiffusion() override {
		DirectionalDiffusion*  phi_diffusion
			= new DirectionalDiffusion(species_indexes_map_["Phi"], species_indexes_map_["Phi"],
				diffusion_coff, bias_diffusion_coff, bias_direction);
		species_diffusion_.push_back(phi_diffusion);
	};
};


/**
 * application dependent solid body initial condition
 */
class ThermosolidBodyInitialCondition
	: public  DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{
		
		if (-BW <= pos_n_[index_i][1] && pos_n_[index_i][1] <= 0.0)
		{
			species_n_[phi_][index_i] = phi_lower_wall;
		}

		if (DH <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH+BW)
		{
			species_n_[phi_][index_i] = phi_upper_wall;
		}
		
	};
public: 
	ThermosolidBodyInitialCondition(SolidBody* diffusion_solid_body)
		: DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>(diffusion_solid_body) {
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};

/**
 * application dependent fluid body initial condition
 */
class ThermofluidBodyInitialCondition
	: public  DiffusionReactionInitialCondition< FluidBody, FluidParticles, WeaklyCompressibleFluid>
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
	ThermofluidBodyInitialCondition(FluidBody* diffusion_fluid_body)
		: DiffusionReactionInitialCondition<FluidBody, FluidParticles, WeaklyCompressibleFluid >(diffusion_fluid_body) {
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};

/**
 *Set thermal relaxation between different bodies 
 */
class ThermalRelaxationComplex
	: public RelaxationOfAllDiffusionSpeciesRK2<FluidBody, FluidParticles, WeaklyCompressibleFluid,
	RelaxationOfAllDiffussionSpeciesComplex<FluidBody, FluidParticles, WeaklyCompressibleFluid, SolidBody, SolidParticles, Solid>,
	SPHBodyComplexRelation>
{
public:
	ThermalRelaxationComplex(SPHBodyComplexRelation* body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation) {};
	virtual ~ThermalRelaxationComplex() {};
};

/** Case dependent inflow boundary condition. */
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
		t_ref = 2.0;
	}
	Vecd getTargetVelocity(Vecd& position, Vecd& velocity)
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0) {
			u = 6.0* u_ave_ * position[1] * (DH - position[1]) / DH / DH;
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

/** an observer body to measure temperature at given poistions */
class TemperatureObserver : public FictitiousBody
{
public:
	TemperatureObserver(SPHSystem& system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
	{
		/** A measuring point at the center of the channel */
			Vec2d point_coordinate(0.0,DH*0.5);
			body_input_points_volumes_.push_back(make_pair(point_coordinate, 0.0));
			
	}
};


/** The main program. */
int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref, 4);
	GlobalStaticVariables::physical_time_ = 0.0;

	/**
	 * @brief Creating body, materials and particles for a ThermofluidBody .
	 */
	ThermofluidBody *thermofluid_body = new ThermofluidBody(system, "ThermofluidBody", 0);
	ThermofluidBodyMaterial *thermofluid_body_material = new ThermofluidBodyMaterial();
	DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>	
		diffusion_fluid_body_particles(thermofluid_body, thermofluid_body_material);

	/**
    * @brief Creating body and particles for the ThermosolidBody.
    */
	ThermosolidBody *thermosolid_body = new ThermosolidBody(system, "ThermosolidBody", 0);
	ThermosolidBodyMaterial *thermosolid_body_material = new ThermosolidBodyMaterial();
	DiffusionReactionParticles<SolidParticles, Solid>	
		diffusion_solid_body_particles(thermosolid_body, thermosolid_body_material);

	/**
	 * @brief 	Particle and body creation of fluid observers.
	 */
	TemperatureObserver* temperature_observer = new TemperatureObserver(system, "FluidObserver", 0);
	BaseParticles	temperature_observer_particles(temperature_observer);
	
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToVtu 				write_real_body_states(in_output, system.real_bodies_);

	/** topology */
	SPHBodyInnerRelation* fluid_body_inner = new SPHBodyInnerRelation(thermofluid_body);
	SPHBodyInnerRelation* solid_body_inner = new SPHBodyInnerRelation(thermosolid_body);
	SPHBodyComplexRelation* fluid_body_complex = new SPHBodyComplexRelation(fluid_body_inner, {thermosolid_body });
	SPHBodyContactRelation* fluid_observer_contact = new SPHBodyContactRelation(temperature_observer, {thermofluid_body});
	
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition(thermofluid_body, 0);

	/**
	 * The main dynamics algorithm is defined start here.
	 */

	 /** Case setup */
	ThermosolidBodyInitialCondition thermosolid_condition(thermosolid_body);
	ThermofluidBodyInitialCondition thermofluid_initial_condition(thermofluid_body);

	/** Corrected strong configuration for diffusion solid body. */
	solid_dynamics::CorrectConfiguration 			correct_configuration(solid_body_inner);
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(thermofluid_body);

	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_density(fluid_body_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step(thermofluid_body, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step(thermofluid_body);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<FluidBody, FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step(thermofluid_body);
	/** Diffusion process between two diffusion bodies. */
	ThermalRelaxationComplex 	thermal_relaxation_complex(fluid_body_complex);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationFirstHalf
		pressure_relaxation_first_half(fluid_body_complex);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(fluid_body_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAcceleration 	viscous_acceleration(fluid_body_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityFormulation 	transport_velocity_formulation(fluid_body_complex);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInFluidField 	compute_vorticity(fluid_body_inner);
	/** Inflow boundary condition. */
	ParabolicInflow	parabolic_inflow(thermofluid_body, new InflowBuffer(thermofluid_body, "Buffer"));

	/**
	 * @brief Write observation data into files.
	 */
	WriteObservedDiffusionReactionQuantity<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>>
		write_fluid_phi("Phi", in_output, fluid_observer_contact);

	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::vel_n_>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);

	 /** Pre-simultion*/
	 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build up
      * but before the configuration build up. */
	periodic_condition.update_cell_linked_list_.parallel_exec();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the wall. */
	diffusion_solid_body_particles.initializeNormalDirectionFromGeometry();
	correct_configuration.parallel_exec();
	thermosolid_condition.parallel_exec();
	thermofluid_initial_condition.parallel_exec();
	Real dt_thmermal = get_thermal_time_step.parallel_exec();

	/** Output global basic parameters. */
    write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	
	Real End_Time = 10;
	Real dt = 0.0;
	Real D_Time = End_Time / 100.0;	/**< time stamps for output,WriteToFile*/ 
	Real Dt = 0.0;					/**< Default advection time step sizes for fluid. */
	size_t inner_ite_dt = 0;
	int number_of_iterations = 0;
	int screen_output_interval = 40;

	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	/** Main loop starts here. */ 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
 		while (integration_time < D_Time) {
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step.parallel_exec();
			update_fluid_density.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);

			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				dt = SMIN(SMIN(dt_thmermal, get_fluid_time_step.parallel_exec()), Dt - relaxation_time);
				pressure_relaxation_first_half.parallel_exec(dt);
				pressure_relaxation_second_half.parallel_exec(dt);			
				thermal_relaxation_complex.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.exec();
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_condition.bounding_.parallel_exec();
			thermofluid_body->updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
            fluid_body_complex->updateConfiguration();
		}
		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		fluid_observer_contact->updateConfiguration();
		write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	    write_fluid_phi.WriteToFile(GlobalStaticVariables::physical_time_);
	    write_fluid_velocity.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
