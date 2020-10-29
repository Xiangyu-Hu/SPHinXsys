/**
 * @file 	filling_tank.cpp
 * @brief 	2D example to show that a tank is filled by emitter.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding how emitter inflow is working.
 * @author 	Chi Zhang and Xiangyu Hu
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
Real DL = 5.366; 						/**< Tank length. */
Real DH = 5.366; 						/**< Tank height. */
Real particle_spacing_ref = 0.025; 		/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; 	/**< Extending width for wall BCs. */
Real LL = 2.0*BW; 						/**< Inflow region length. */
Real LH = 0.125; 						/**< Inflows region height. */
Real inlet_height = 1.0;				/**< Inflow location height */
Real inlet_distance = - BW;				/**< Inflow location distance */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;						/**< Reference density of fluid. */
Real gravity_g = 1.0;					/**< Gravity force of fluid. */
Real U_f = 2.0*sqrt(gravity_g* (inlet_height + LH));	/**< Characteristic velocity. */
Real c_f = 10.0*U_f;					/**< Reference sound speed. */

/** @brief create a water block shape for the inlet. */
std::vector<Point> CreatWaterBlockShape()
{
	//geometry
	std::vector<Point> water_block_shape;
	water_block_shape.push_back(Point(inlet_distance, inlet_height));
	water_block_shape.push_back(Point(inlet_distance, LH + inlet_height));
	water_block_shape.push_back(Point(LL + inlet_distance, LH + inlet_height));
	water_block_shape.push_back(Point(LL + inlet_distance, inlet_height));
	water_block_shape.push_back(Point(inlet_distance, inlet_height));

	return water_block_shape;
}

/** @brief 	Fluid body definition. */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, string body_name,	int refinement_level)
		: FluidBody(system, body_name, refinement_level)
	{
		/** Geometry definition. */
		std::vector<Point> water_bock_shape = CreatWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_bock_shape, ShapeBooleanOps::add);

		/** Replace the default kernel functions */
		ReplaceKernelFunction(new KernelTabulated<KernelWendlandC2>(smoothing_length_, 20));
	}
};
/**
 * @brief 	Case dependent material properties definition.
 */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		rho_0_ = rho0_f;
		c_0_ = c_f;

		assignDerivedMaterialParameters();
	}
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geometry definition. */
		std::vector<Point> outer_wall_shape;
		outer_wall_shape.push_back(Point(-BW, -BW));
		outer_wall_shape.push_back(Point(-BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, -BW));
		outer_wall_shape.push_back(Point(-BW, -BW));

		std::vector<Point> inner_wall_shape;
		inner_wall_shape.push_back(Point(0.0, 0.0));
		inner_wall_shape.push_back(Point(0.0, DH));
		inner_wall_shape.push_back(Point(DL, DH));
		inner_wall_shape.push_back(Point(DL, 0.0));
		inner_wall_shape.push_back(Point(0.0, 0.0));
		/** Inlet */
		std::vector<Point> water_block_shape = CreatWaterBlockShape();

		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::sub);
	}
};

/**
* @brief constrain the beam base
*/
class Inlet : public BodyPartByParticle
{
public:
	Inlet(FluidBody* fluid_body, string constrained_region_name)
		: BodyPartByParticle(fluid_body, constrained_region_name)
	{
		/** Geometry definition. */
		std::vector<Point> water_block_shape = CreatWaterBlockShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);

		/**  Tag the constrained particle. */
		tagBodyPart();
	}
};

/** inlet inflow condition. */
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

/**
 * @brief 	Fluid observer body definition.
 */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem &system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
	{
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.2), 0.0));
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
	SPHSystem system(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", 0);
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall",	0);
	SolidParticles	wall_particles(wall_boundary);
	/**
	 * @brief 	Particle and body creation of fluid observer.
	 */
	FluidObserver *fluid_observer = new FluidObserver(system, "Fluidobserver", 0);
	BaseParticles 	observer_particles(fluid_observer);
	/** topology */
	SPHBodyComplexRelation* water_block_complex_relation = new SPHBodyComplexRelation(water_block, { wall_boundary });
	SPHBodyContactRelation* fluid_observer_contact_relation = new SPHBodyContactRelation(fluid_observer, { water_block });
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Define external force. */
	Gravity 							gravity(Vecd(0.0, -gravity_g));
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeATimeStep 	initialize_a_fluid_step(water_block, &gravity);
	/** Inlet condition. */
	Inlet* inlet = new Inlet(water_block, "Inlet");
	InletInflowCondition inflow_condition(water_block, inlet);
	fluid_dynamics::EmitterInflowInjecting inflow_emitter(water_block, inlet, 300, 0, true);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensityBySummationFreeSurface 		update_fluid_density(water_block_complex_relation);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 			get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation algorithm by using position verlet time stepping. */
	fluid_dynamics::PressureRelaxationFirstHalfRiemann 
		pressure_relaxation_first_half(water_block_complex_relation);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann 
		pressure_relaxation_second_half(water_block_complex_relation);
	/**
	 * @brief Output.
	 */
	In_Output in_output(system);
	/** Output the body states. */
	WriteBodyStatesToVtu 		write_body_states(in_output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	ReadRestart		read_restart_files(in_output, system.real_bodies_);
	WriteRestart	write_restart_files(in_output, system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	WriteTotalMechanicalEnergy 	write_water_mechanical_energy(in_output, water_block, &gravity);
	/** output the observed data from fluid body. */
	WriteAnObservedQuantity<Real, FluidParticles, &FluidParticles::p_>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact_relation);
	
	/**
	 * @brief Setup configurations and initial conditions.
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	/**
	 * @brief The time stepping starts here.
	 */
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = read_restart_files.ReadRestartFiles(system.restart_step_);
		water_block->updateCellLinkedList();
		water_block_complex_relation->updateConfiguration();
	}
	/** Output the start states of bodies. */
	write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Output the Hydrostatic mechanical energy of fluid. */
	write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval*10;
	Real End_Time = 50.0; 	/**< End time. */
	Real D_Time = 0.1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	/**
	 * @brief 	Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();

			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation_first_half.parallel_exec(dt);
				inflow_condition.parallel_exec();
				pressure_relaxation_second_half.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					write_restart_files.WriteToFile(Real(number_of_iterations));
			}
			number_of_iterations++;

			/** impose inflow condition*/
			inflow_emitter.exec();

			/** Update cell linked list and configuration. */
			water_block->updateCellLinkedList();
			water_block_complex_relation->updateConfiguration();
			fluid_observer_contact_relation->updateConfiguration();
		}

		tick_count t2 = tick_count::now();
		write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
		write_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
		write_recorded_water_pressure.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;

	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << endl;

	return 0;
}
