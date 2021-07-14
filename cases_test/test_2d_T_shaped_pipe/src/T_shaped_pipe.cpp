/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"

using namespace SPH;
/**
 * @brief Basic geometry parameters.
 */
Real DL = 5.0; 						    /**< Domain length. */
Real DH = 3.0; 						    /**< Domain height. */
Real resolution_ref = 0.15; 			/**< Initial reference particle spacing. */
Real BW = resolution_ref * 4; 			/**< Sponge region to impose injection. */
Real DL_sponge = resolution_ref * 20; 	/**< Sponge region to impose emitter. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -DH), Vec2d(DL + BW, 2.0 * DH));
/** Prescribed fluid body domain bounds*/
BoundingBox fluid_body_domain_bounds(Vec2d(-DL_sponge, -DH), Vec2d(DL + BW, 2.0 * DH));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real U_f = 1.0;	                        /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;					/**< Reference sound speed. */
Real Re = 100.0;		                /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re;	    /**< Dynamics viscosity. */
/**
* @brief define fluid
*/
/** create a water block shape. */
std::vector<Vecd> CreateWaterBlockShape()
{
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, DH));
	water_block_shape.push_back(Vecd(0.7 * DL, DH));
	water_block_shape.push_back(Vecd(0.7 * DL, 2.0 * DH));
	water_block_shape.push_back(Vecd(DL, 2.0 * DH));
	water_block_shape.push_back(Vecd(DL, -DH));
	water_block_shape.push_back(Vecd(0.7 * DL, -DH));
	water_block_shape.push_back(Vecd(0.7 * DL, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

	return water_block_shape;
}
/** @brief 	fluid body definition. */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& system, string body_name)
		: FluidBody(system, body_name)
	{
		std::vector<Vecd> water_bock_shape = CreateWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_bock_shape, ShapeBooleanOps::add);
	}
};
/** Case fluid material properties definition. */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		rho0_ = rho0_f;
		c0_ = c_f;
		mu_ = mu_f;

		assignDerivedMaterialParameters();
	}
};
/**
* @brief define wall boundary
*/
/** create a outer wall polygen. */
std::vector<Vecd> CreateOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-DL_sponge, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge, DH + BW));
	outer_wall_shape.push_back(Vecd(0.7 * DL - BW, DH + BW));
	outer_wall_shape.push_back(Vecd(0.7 * DL - BW, 2.0 * DH));
	outer_wall_shape.push_back(Vecd(DL + BW, 2.0 * DH));
	outer_wall_shape.push_back(Vecd(DL + BW, -DH));
	outer_wall_shape.push_back(Vecd(0.7 * DL - BW, -DH));
	outer_wall_shape.push_back(Vecd(0.7 * DL - BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge, -BW));

	return outer_wall_shape;
}
/** create a inner wall polygen. */
std::vector<Vecd> CreateInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, DH));
	inner_wall_shape.push_back(Vecd(0.7 * DL, DH));
	inner_wall_shape.push_back(Vecd(0.7 * DL, 2.0 * DH + BW));
	inner_wall_shape.push_back(Vecd(DL, 2.0 * DH + BW));
	inner_wall_shape.push_back(Vecd(DL, -DH - BW));
	inner_wall_shape.push_back(Vecd(0.7 * DL, -DH - BW));
	inner_wall_shape.push_back(Vecd(0.7 * DL, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - BW, 0.0));

	return inner_wall_shape;
}
/** Wall boundary body definition. */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		std::vector<Vecd> outer_wall_shape = CreateOuterWallShape();
		std::vector<Vecd> inner_wall_shape = CreateInnerWallShape();
		std::vector<Vecd> water_block_shape = CreateWaterBlockShape();

		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	}
};
/**
* @brief define the emitter buffer
*/
/** create the emitter buffer shape . */
std::vector<Vecd> CreateEmitterBufferShape()
{
	std::vector<Vecd> emitter_buffer_shape;
	emitter_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
	emitter_buffer_shape.push_back(Vecd(-DL_sponge, DH));
	emitter_buffer_shape.push_back(Vecd(0.0, DH));
	emitter_buffer_shape.push_back(Vecd(0.0, 0.0));
	emitter_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

	return emitter_buffer_shape;
}
/** define emitter buffer body*/
class EmitterBuffer : public BodyPartByParticle
{
public:
	EmitterBuffer(FluidBody* fluid_body, std::string constrained_region_name)
		: BodyPartByParticle(fluid_body, constrained_region_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> emitter_buffer_shape = CreateEmitterBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(emitter_buffer_shape, ShapeBooleanOps::add);

		/**  Tag the constrained particle. */
		tagBodyPart();
	}
};
/** define emitter inflow boundary conditon*/
class EmitterInflowCondition : public fluid_dynamics::InletOutletInflowCondition
{
	Real u_ave_, u_ref_, t_ref_;
public:
	EmitterInflowCondition(FluidBody* body, BodyPartByParticle* body_part)
		: InletOutletInflowCondition(body, body_part)
	{
		u_ave_ = 0.0;
		u_ref_ = U_f;
		t_ref_ = 2.0;
		SetInflowParameters();
	}

	Vecd getTargetVelocity(Vecd& position, Vecd& velocity) override
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
		u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
	}

	void SetInflowParameters() override
	{
		inflow_pressure_ = 0.0;
	}
};
/**
* @brief define the injection buffer
*/
/** create the injection buffer shape. */
std::vector<Vecd> CreatInjectionBufferShape()
{
	std::vector<Vecd> injection_buffer_shape;
	injection_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
	injection_buffer_shape.push_back(Vecd(-DL_sponge, DH));
	injection_buffer_shape.push_back(Vecd(-DL_sponge + BW, DH));
	injection_buffer_shape.push_back(Vecd(-DL_sponge + BW, 0.0));
	injection_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

	return injection_buffer_shape;
}
/** define injection buffer */
class InjectionBuffer : public BodyPartByParticle
{
public:
	InjectionBuffer(FluidBody* fluid_body, string constrained_region_name)
		: BodyPartByParticle(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> emitter_buffer_shape = CreatInjectionBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(emitter_buffer_shape, ShapeBooleanOps::add);

		//tag the constrained particle
		tagBodyPart();
	}
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char* av[])
{
	/** Build up a SPHSystem */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);
	In_Output in_output(system);

	/**
	 * @brief Creating body, materials and particles for a water block.
	 */
	WaterBlock* water_block = new WaterBlock(system, "WaterBody");
	water_block->setBodyDomainBounds(fluid_body_domain_bounds);
	WaterMaterial* water_material = new WaterMaterial();
	FluidParticles fluid_particles(water_block, water_material);
	/**
	 * @brief Creating body, materials and particles for wall boundary.
	 */
	WallBoundary* wall_boundary = new WallBoundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	/**
    * @brief File Output.
    */
	BodyStatesRecordingToVtu write_body_states(in_output, system.real_bodies_);
	RestartIO restart_io(in_output, system.real_bodies_);
	/**
	 * @brief Define body relation map.
	 * @brief The contact map gives the topological connections between the bodies.
	 * @brief Basically the the range of bodies to build neighbor particle lists.
	 */
	BodyRelationInner* water_block_inner = new BodyRelationInner(water_block);
	ComplexBodyRelation* water_block_complex_relation = new ComplexBodyRelation(water_block_inner, { wall_boundary });
	/**
	 * @brief Define all numerical methods which are used in this case.
	 */
	 /** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block);
	/** Emitter condition. */
	EmitterBuffer* emitter = new EmitterBuffer(water_block, "Emitter");
	EmitterInflowCondition inflow_condition(water_block, emitter);
	/** Injection condition. */
	InjectionBuffer* injection = new InjectionBuffer(water_block, "Injection");
	fluid_dynamics::EmitterInflowInjecting inflow_emitter(water_block, injection, 300, 0, true);
	/** time-space method to detect surface particles. */
	fluid_dynamics::SurfaceParticlesIndicator
		inlet_outlet_surface_particle_indicator(water_block_complex_relation, water_block_inner);
	/** Evaluation of density by freestream approach. */
	fluid_dynamics::DensitySummationFreeStreamComplex update_density_by_summation(water_block_complex_relation);
	/** We can output a method-specific particle data for debug reason */
	fluid_particles.addAVariableToWrite<indexScalar, Real>("Pressure");
	fluid_particles.addAVariableToWrite<indexInteger, int>("SurfaceIndicator");
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex_relation);
	/** Density relaxation. */
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex_relation);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex_relation);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex_relation);
	/** recycle real fluid particle to buffer particles at outlet. */
	OpenBoundaryConditionInAxisDirection tansfer_to_buffer_particles_lower_bound(water_block, yAxis, negativeDirection);
	OpenBoundaryConditionInAxisDirection tansfer_to_buffer_particles_upper_bound(water_block, yAxis, positiveDirection);
	/**
	* @brief Setup computing and initial conditions.
	*/
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	inlet_outlet_surface_particle_indicator.first_layer.parallel_exec();
	inlet_outlet_surface_particle_indicator.second_and_third_layers.parallel_exec();
	/**
	* @brief The time stepping starts here.
	*/
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		water_block->updateCellLinkedList();
		water_block_complex_relation->updateConfiguration();
	}
	write_body_states.writeToFile(GlobalStaticVariables::physical_time_);
	/**
	* @brief Time stepping control parameters.
	*/
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 100.0; 	/**< End time. */
	Real D_Time = End_Time / 200.0; /**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				pressure_relaxation.parallel_exec(dt);
				inflow_condition.parallel_exec();
				density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** impose inflow condition*/
			inflow_emitter.exec();
			tansfer_to_buffer_particles_lower_bound.particle_type_transfer.parallel_exec();
			tansfer_to_buffer_particles_upper_bound.particle_type_transfer.parallel_exec();

			/** Update cell linked list and configuration. */
			water_block->updateCellLinkedList();
			water_block_complex_relation->updateConfiguration();
			inlet_outlet_surface_particle_indicator.first_layer.parallel_exec();
			inlet_outlet_surface_particle_indicator.second_and_third_layers.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		write_body_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;

	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;

	return 0;
}
