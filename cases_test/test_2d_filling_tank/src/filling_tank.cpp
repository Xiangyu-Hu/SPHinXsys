/**
 * @file 	filling_tank.cpp
 * @brief 	2D example to show that a tank is filled by emitter.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding how emitter inflow is working.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
 //----------------------------------------------------------------------
 //	Global geometry, material parameters and numerical setup.
 //----------------------------------------------------------------------
Real DL = 5.366; 						/**< Tank length. */
Real DH = 5.366; 						/**< Tank height. */
Real resolution_ref = 0.025; 			/**< Initial reference particle spacing. */
Real BW = resolution_ref * 4; 			/**< Extending width for wall boundarys. */
Real LL = 2.0*BW; 						/**< Inflow region length. */
Real LH = 0.125; 						/**< Inflows region height. */
Real inlet_height = 1.0;				/**< Inflow location height */
Real inlet_distance = - BW;				/**< Inflow location distance */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real gravity_g = 1.0;					/**< Gravity force of fluid. */
Real U_f = 2.0*sqrt(gravity_g* (inlet_height + LH));	/**< Characteristic velocity. */
Real c_f = 10.0*U_f;					/**< Reference sound speed. */
/** create a outer wall polygen. */
std::vector<Vecd> CreateOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));

	return outer_wall_shape;
}
/** create a inner wall polygen. */
std::vector<Vecd> CreateInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

	return inner_wall_shape;
}
/** create a water block shape for the inlet. */
std::vector<Vecd> CreateWaterBlockShape()
{
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(inlet_distance, inlet_height));
	water_block_shape.push_back(Vecd(inlet_distance, LH + inlet_height));
	water_block_shape.push_back(Vecd(LL + inlet_distance, LH + inlet_height));
	water_block_shape.push_back(Vecd(LL + inlet_distance, inlet_height));
	water_block_shape.push_back(Vecd(inlet_distance, inlet_height));

	return water_block_shape;
}
/** @brief 	Fluid body definition. */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, std::string body_name, ParticleAdaptation* particle_adaptation)
		: FluidBody(system, body_name, particle_adaptation)
	{
		/** initial water block is the inlet */
		std::vector<Vecd> water_bock_shape = CreateWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_bock_shape, ShapeBooleanOps::add);
	}
};
/** Case dependent material properties definition. */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		rho0_ = rho0_f;
		c0_ = c_f;

		assignDerivedMaterialParameters();
	}
};
/** Wall boundary body definition. */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape = CreateOuterWallShape();
		std::vector<Vecd> inner_wall_shape = CreateInnerWallShape();
		std::vector<Vecd> water_block_shape = CreateWaterBlockShape();

		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	}
};
/** Inelt as a bodypart by particle */
class Inlet : public BodyPartByParticle
{
public:
	Inlet(FluidBody* fluid_body, std::string constrained_region_name)
		: BodyPartByParticle(fluid_body, constrained_region_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_block_shape = CreateWaterBlockShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);

		/**  Tag the constrained particle. */
		tagBodyPart();
	}
};
/** Inlet inflow condition. */
class InletInflowCondition : public fluid_dynamics::EmitterInflowCondition
{
public:
	InletInflowCondition(FluidBody* body, BodyPartByParticle* body_part)
		: EmitterInflowCondition(body, body_part)
	{
		SetInflowParameters();
	}
	
	Vecd getTargetVelocity(Vecd& position, Vecd& velocity) override {
		return Vec2d(2.0, 0.0);
	}
	void SetInflowParameters() override {
		inflow_pressure_ = 0.0;
	}
};
/** Fluid observer body definition. */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem &system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		body_input_points_volumes_.push_back(std::make_pair(Vecd(DL, 0.2), 0.0));
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	/** Build up a SPHSystem */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	ParticleAdaptation* particle_adaptation = new ParticleAdaptation();
	particle_adaptation->replaceKernel(new KernelTabulated<KernelWendlandC2>(20));
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", particle_adaptation);
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**note that, as particle sort is activated (by default) for fluid particles, 
	 * the output occasionally does not reflect the real free surface indication due to sorting. */
	fluid_particles.addAVariableToWrite<indexInteger, int>("SurfaceIndicator");

	WallBoundary *wall_boundary = new WallBoundary(system, "Wall");
	SolidParticles	wall_particles(wall_boundary);

	FluidObserver *fluid_observer = new FluidObserver(system, "Fluidobserver");
	BaseParticles 	observer_particles(fluid_observer);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexBodyRelation* water_block_complex = new ComplexBodyRelation(water_block, { wall_boundary });
	BodyRelationContact* fluid_observer_contact_relation = new BodyRelationContact(fluid_observer, { water_block });
	//----------------------------------------------------------------------
	//	Define all numerical methods which are used in this case.
	//----------------------------------------------------------------------
	Gravity 				gravity(Vecd(0.0, -gravity_g));
	TimeStepInitialization 	initialize_a_fluid_step(water_block, &gravity);
	Inlet* inlet = new Inlet(water_block, "Inlet");
	InletInflowCondition inflow_condition(water_block, inlet);
	fluid_dynamics::EmitterInflowInjecting inflow_emitter(water_block, inlet, 300, 0, true);
	fluid_dynamics::DensitySummationFreeSurfaceComplex	update_density_by_summation(water_block_complex);
	fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex	indicate_free_surface(water_block_complex);
	/** We can output a method-specific particle data for debug reason */
	fluid_particles.addAVariableToWrite<indexScalar, Real>("PositionDivergence");
	fluid_particles.addAVariableToWrite<indexInteger, int>("SurfaceIndicator");
	fluid_dynamics::AdvectionTimeStepSize 			get_fluid_advection_time_step_size(water_block, U_f);
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	//----------------------------------------------------------------------
	//	File Output
	//----------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtu 		body_states_recording(in_output, system.real_bodies_);
	RestartIO		restart_io(in_output, system.real_bodies_);
	BodyReducedQuantityRecording<TotalMechanicalEnergy> 	
		write_water_mechanical_energy(in_output, water_block, &gravity);
	ObservedQuantityRecording<indexScalar, Real>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact_relation);
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	indicate_free_surface.parallel_exec();
	//----------------------------------------------------------------------
	//	The time stepping starts here.
	//----------------------------------------------------------------------
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		water_block->updateCellLinkedList();
		water_block_complex->updateConfiguration();
	}
	body_states_recording.writeToFile(0);
	write_water_mechanical_energy.writeToFile(0);
	//----------------------------------------------------------------------
	//	Time stepping control parameters.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval*10;
	Real End_Time = 30.0; 	/**< End time. */
	Real D_Time = 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();			
			update_density_by_summation.parallel_exec();

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation.parallel_exec(dt);
				inflow_condition.parallel_exec();
				density_relaxation.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** impose inflow condition*/
			inflow_emitter.exec();

			/** Update cell linked list and configuration. */
			water_block->updateCellLinkedList();
			water_block_complex->updateConfiguration();
			fluid_observer_contact_relation->updateConfiguration();
			indicate_free_surface.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		write_water_mechanical_energy.writeToFile(number_of_iterations);
		body_states_recording.writeToFile();
		write_recorded_water_pressure.writeToFile(number_of_iterations);
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
