/**
 * @file 	flow_wall_heat_transfer.cpp
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



/**  Diffusion fluid body definition */
class DiffusionFluidBody : public FluidBody
{
public: 
	DiffusionFluidBody(SPHSystem &system, string body_name, int refinement_level)
		: FluidBody(system, body_name, refinement_level)
	{	
		std::vector<Point> body_shape = CreatShape();	
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(body_shape, ShapeBooleanOps::add);
	}
};

/**  Diffusion solid body definition */
class DiffusionSolidBody : public SolidBody
{
public:
	DiffusionSolidBody(SPHSystem &system, string body_name, int refinement_level)
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
 * Setup diffusion material properties for diffusion fluid body 
 */
class DiffusionFluidBodyMaterial
	: public DiffusionReactionMaterial<FluidParticles, WeaklyCompressibleFluid>
{
public:
	DiffusionFluidBodyMaterial()
		: DiffusionReactionMaterial<FluidParticles, WeaklyCompressibleFluid>()
	{
		rho_0_ = rho0_f;
		c_0_ = c_f;
		mu_ = mu_f;

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
 * Setup diffusion material properties for diffusion solid body
 */
class DiffusionSolidBodyMaterial
	: public DiffusionReactionMaterial<SolidParticles, Solid>
{
public:
	DiffusionSolidBodyMaterial()
		: DiffusionReactionMaterial<SolidParticles, Solid>()
	{
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
class DiffusionSolidBodyInitialCondition
	: public  DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{
		
		if (-BW <= pos_n_[index_i][1] && pos_n_[index_i][1] <= 0.0)
		{
			species_n_[phi_][index_i] = 40.0;
		}

		if (DH <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH+BW)
		{
			species_n_[phi_][index_i] = 20.0;
		}
		
	};
public: 
	DiffusionSolidBodyInitialCondition(SolidBody* diffusion_solid_body)
		: DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>(diffusion_solid_body) {
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};

/**
 * application dependent solid body initial condition
 */
class DiffusionFluidBodyInitialCondition
	: public  DiffusionReactionInitialCondition< FluidBody, FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{

		if (0 <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH)
		{
			species_n_[phi_][index_i] = 20.0;
		}

	};
public:
	DiffusionFluidBodyInitialCondition(FluidBody* diffusion_fluid_body)
		: DiffusionReactionInitialCondition<FluidBody, FluidParticles, WeaklyCompressibleFluid >(diffusion_fluid_body) {
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};

/**
 *Set diffusion ComplexRelaxation between different bodies 
 */
class DiffusionBodyComplexRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<FluidBody, FluidParticles, WeaklyCompressibleFluid,
	RelaxationOfAllDiffussionSpeciesComplex<FluidBody, FluidParticles, WeaklyCompressibleFluid, SolidBody, SolidParticles, Solid>,
	SPHBodyComplexRelation>
{
public:
	DiffusionBodyComplexRelaxation(SPHBodyComplexRelation* body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation) {};
	virtual ~DiffusionBodyComplexRelaxation() {};
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

/** an observer body to measure phi dependent of temperature */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem& system, string body_name, int refinement_level)
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
	 * @brief Creating body, materials and particles for a DiffusionFluidBody .
	 */
	DiffusionFluidBody *diffusion_fluid_body = new DiffusionFluidBody(system, "DiffusionFluidBody", 0); 
	DiffusionFluidBodyMaterial *diffusion_fluid_body_material = new DiffusionFluidBodyMaterial();
	DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>	diffusion_fluid_body_particles(diffusion_fluid_body, diffusion_fluid_body_material);

	/**
    * @brief Creating body and particles for the DiffusionSolidBody.
    */
	DiffusionSolidBody *diffusion_solid_body = new DiffusionSolidBody(system, "DiffusionSolidBody", 0);
	DiffusionSolidBodyMaterial *diffusion_solid_body_material = new DiffusionSolidBodyMaterial();
	DiffusionReactionParticles<SolidParticles, Solid>	diffusion_solid_body_particles(diffusion_solid_body, diffusion_solid_body_material);

	/**
	 * @brief 	Particle and body creation of fluid observers.
	 */
	FluidObserver* fluid_observer = new FluidObserver(system, "FluidObserver", 0);
	BaseParticles				flow_observer_particles(fluid_observer);
	
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToVtu 				write_real_body_states(in_output, system.real_bodies_);

	/** topology */
	SPHBodyInnerRelation* diffusion_fluid_body_inner_relation = new SPHBodyInnerRelation(diffusion_fluid_body);
	SPHBodyInnerRelation* diffusion_solid_body_inner_relation = new SPHBodyInnerRelation(diffusion_solid_body);
	SPHBodyComplexRelation* diffusion_fluid_body_complex = new SPHBodyComplexRelation(diffusion_fluid_body_inner_relation, {diffusion_solid_body });
	SPHBodyContactRelation* fluid_observer_contact = new SPHBodyContactRelation(fluid_observer, {diffusion_fluid_body});
	
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition(diffusion_fluid_body, 0);

	/**
	 * The main dynamics algorithm is defined start here.
	 */

	 /** Case setup */
	DiffusionSolidBodyInitialCondition setup_diffusion_initial_condition(diffusion_solid_body);
	DiffusionFluidBodyInitialCondition setup_fluid_diffusion_initial_condition(diffusion_fluid_body);

	/** Corrected strong configuration for diffusion solid body. */
	solid_dynamics::CorrectConfiguration 			correct_configuration(diffusion_solid_body_inner_relation);
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(diffusion_fluid_body);

	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummation 			update_fluid_density(diffusion_fluid_body_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step_size(diffusion_fluid_body, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(diffusion_fluid_body);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationFirstHalf
		pressure_relaxation_first_half(diffusion_fluid_body_complex);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(diffusion_fluid_body_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAcceleration 	viscous_acceleration(diffusion_fluid_body_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityFormulation 	transport_velocity_formulation(diffusion_fluid_body_complex);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInFluidField 	compute_vorticity(diffusion_fluid_body_inner_relation);
	/** Inflow boundary condition. */
	ParabolicInflow	parabolic_inflow(diffusion_fluid_body, new InflowBuffer(diffusion_fluid_body, "Buffer"));

	/**
	 * @brief Write observation data into files.
	 */
	WriteObservedDiffusionReactionQuantity<DiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>>
		write_fluid_phi("Phi", in_output, fluid_observer_contact);

	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::vel_n_>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);
	

	/** Diffusion process between two diffusion bodies. */
	DiffusionBodyComplexRelaxation 			diffusion_complex_relaxation(diffusion_fluid_body_complex);
	


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
	setup_diffusion_initial_condition.parallel_exec();
	setup_fluid_diffusion_initial_condition.parallel_exec();

	/** Output global basic parameters. */
    write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	
	Real End_Time = 10;
	Real dt = 0.0;
	Real D_Time = End_Time / 200.0;	/**< time stamps for output,WriteToFile*/ 
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
			
			setup_diffusion_initial_condition.exec();
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);

			inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				/** Fluid pressure relaxation, first half. */
				pressure_relaxation_first_half.parallel_exec(dt);
				/** Fluid pressure relaxation, second half. */
				pressure_relaxation_second_half.parallel_exec(dt);			

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.exec();
				inner_ite_dt++;

				diffusion_complex_relaxation.parallel_exec(dt);

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
			diffusion_fluid_body->updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
            diffusion_fluid_body_complex->updateConfiguration();
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


