/**
 * @file 	heat_transfer.cpp
 * @brief 	This is a test to validate heat transfer between a flow and channel walls.
 * @author 	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 */
/** SPHinXsys Library. */
#include "sphinxsys.h"
using namespace SPH;


/** Geometry parameter. */
Real DL = 2.0;                /**< Channel length. */
Real DH = 0.4;              /**< Channel height. */
Real resolution_ref = DH / 25.0;          /**< Global reference reoslution. */
Real DL_sponge = resolution_ref * 20.0;	/**< Sponge region to impose inflow condition. */
/**< Boundary width, determined by specific layer of boundary particles. */
Real BW = resolution_ref * 4.0; 	/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));

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


/** create a water block buffer shape. */
std::vector<Vecd> CreatInflowBufferShape()
{
	std::vector<Vecd> inflow_buffer_shape;
	inflow_buffer_shape.push_back(Vecd(0.0 - DL_sponge, 0.0));
	inflow_buffer_shape.push_back(Vecd(0.0 - DL_sponge, DH));
	inflow_buffer_shape.push_back(Vecd(0, DH));
	inflow_buffer_shape.push_back(Vecd(0, 0));
	inflow_buffer_shape.push_back(Vecd(0.0 - DL_sponge, 0.0));

	return inflow_buffer_shape;
}


/** inflow buffer */
class InflowBuffer : public BodyPartByCell
{
public:
	InflowBuffer(FluidBody* fluid_body, std::string constrained_region_name)
		: BodyPartByCell(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> inflow_buffer_shape = CreatInflowBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
		tagBodyPart();
	}
};



/**  Thermofluid body definition */
class ThermofluidBody : public FluidBody
{
public: 
	ThermofluidBody(SPHSystem &system, std::string body_name)
		: FluidBody(system, body_name)
	{	
		std::vector<Vecd> body_shape = createShape();	
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(body_shape, ShapeBooleanOps::add);
	}
};

/**  Thermosolid body definition */
class ThermosolidBody : public SolidBody
{
public:
	ThermosolidBody(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		std::vector<Vecd>  outer_wall_shape = createOuterWallShape();
		std::vector<Vecd> inner_wall_shape = createInnerWallShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);

	}
};


/**
 * Setup heat conduction material properties for diffusion fluid body 
 */
class ThermofluidBodyMaterial
	: public DiffusionReactionMaterial<FluidParticles, WeaklyCompressibleFluid>
{
public:
	ThermofluidBodyMaterial()
		: DiffusionReactionMaterial<FluidParticles, WeaklyCompressibleFluid>()
	{
		rho0_ = rho0_f;
		c0_ = c_f;
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
	ComplexBodyRelation>
{
public:
	ThermalRelaxationComplex(ComplexBodyRelation* body_complex_relation)
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
	TemperatureObserver(SPHSystem& system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		/** A measuring point at the center of the channel */
			Vec2d point_coordinate(0.0,DH*0.5);
			body_input_points_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
			
	}
};


/** The main program. */
int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(system_domain_bounds, resolution_ref);
	GlobalStaticVariables::physical_time_ = 0.0;

	/**
	 * @brief Creating body, materials and particles for a ThermofluidBody .
	 */
	ThermofluidBody *thermofluid_body = new ThermofluidBody(system, "ThermofluidBody");
	ThermofluidBodyMaterial *thermofluid_body_material = new ThermofluidBodyMaterial();
	DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>	
		diffusion_fluid_body_particles(thermofluid_body, thermofluid_body_material);

	/**
    * @brief Creating body and particles for the ThermosolidBody.
    */
	ThermosolidBody *thermosolid_body = new ThermosolidBody(system, "ThermosolidBody");
	ThermosolidBodyMaterial *thermosolid_body_material = new ThermosolidBodyMaterial();
	DiffusionReactionParticles<SolidParticles, Solid>	
		diffusion_solid_body_particles(thermosolid_body, thermosolid_body_material);

	/**
	 * @brief 	Particle and body creation of fluid observers.
	 */
	TemperatureObserver* temperature_observer = new TemperatureObserver(system, "FluidObserver");
	BaseParticles	temperature_observer_particles(temperature_observer);
	
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	In_Output 							in_output(system);
	BodyStatesRecordingToVtu 				write_real_body_states(in_output, system.real_bodies_);

	/** topology */
	BodyRelationInner* fluid_body_inner = new BodyRelationInner(thermofluid_body);
	BodyRelationInner* solid_body_inner = new BodyRelationInner(thermosolid_body);
	ComplexBodyRelation* fluid_body_complex = new ComplexBodyRelation(fluid_body_inner, {thermosolid_body });
	BodyRelationContact* fluid_observer_contact = new BodyRelationContact(temperature_observer, {thermofluid_body});
	
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition(thermofluid_body, xAxis);

	/**
	 * The main dynamics algorithm is defined start here.
	 */

	 /** Case setup */
	ThermosolidBodyInitialCondition thermosolid_condition(thermosolid_body);
	ThermofluidBodyInitialCondition thermofluid_initial_condition(thermofluid_body);

	/** Corrected strong configuration for diffusion solid body. */
	solid_dynamics::CorrectConfiguration 			correct_configuration(solid_body_inner);
	 /** Initialize particle acceleration. */
	TimeStepInitialization 	initialize_a_fluid_step(thermofluid_body);

	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationComplex	update_density_by_summation(fluid_body_complex);
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
	fluid_dynamics::PressureRelaxationWithWall	pressure_relaxation(fluid_body_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(fluid_body_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall 	viscous_acceleration(fluid_body_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex	transport_velocity_correction(fluid_body_complex);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInner 	compute_vorticity(fluid_body_inner);
	/** Inflow boundary condition. */
	ParabolicInflow	parabolic_inflow(thermofluid_body, new InflowBuffer(thermofluid_body, "Buffer"));

	/**
	 * @brief Write observation data into files.
	 */
	ObservedQuantityRecording<indexScalar, Real>
		write_fluid_phi("Phi", in_output, fluid_observer_contact);
	ObservedQuantityRecording<indexVector, Vecd>
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
	Real dt_thermal = get_thermal_time_step.parallel_exec();

	/** Output global basic parameters. */
    write_real_body_states.writeToFile(0);
	
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
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				dt = SMIN(SMIN(dt_thermal, get_fluid_time_step.parallel_exec()), Dt);
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
			thermofluid_body->updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
            fluid_body_complex->updateConfiguration();
		}
		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		fluid_observer_contact->updateConfiguration();
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

	return 0;
}
