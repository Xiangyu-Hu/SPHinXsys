/**
 * @file 	elastic_gate.cpp
 * @brief 	2D elastic gate deformation due to dam break force.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid-structure-interaction (FSI) simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 500.0;						/**< Tank length. */
Real DH = 200.1;						/**< Tank height. */
Real Dam_L = 100.0;						/**< Water block width. */
Real Dam_H = 140.0;						/**< Water block height. */
Real Gate_width = 5.0;					/**< Width of the gate. */
Real Base_bottom_position = 79.0;		/**< Position of gate base. (In Y direction) */
Real resolution_ref = Gate_width / 2.0; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4.0;			/**< Extending width for BCs. */
/** The offset that the rubber gate shifted above the tank. */
Real dp_s = 0.5 * resolution_ref;
Vec2d offset = Vec2d(0.0, Base_bottom_position - floor(Base_bottom_position / dp_s) * dp_s);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/** Define the corner points of the water block geometry. */
Vec2d DamP_lb(DL - Dam_L, 0.0);	  /**< Left bottom. */
Vec2d DamP_lt(DL - Dam_L, Dam_H); /**< Left top. */
Vec2d DamP_rt(DL, Dam_H);		  /**< Right top. */
Vec2d DamP_rb(DL, 0.0);			  /**< Right bottom. */
/** Define the corner points of the gate geometry. */
Vec2d GateP_lb(DL - Dam_L - Gate_width, 0.0);		 /**< Left bottom. */
Vec2d GateP_lt(DL - Dam_L - Gate_width, Dam_H + BW); /**< Left top. */
Vec2d GateP_rt(DL - Dam_L, Dam_H + BW);				 /**< Right top. */
Vec2d GateP_rb(DL - Dam_L, 0.0);					 /**< Right bottom. */
/** Define the corner points of the gate constrain. */
Vec2d ConstrainP_lb(DL - Dam_L - Gate_width, Base_bottom_position); /**< Left bottom. */
Vec2d ConstrainP_lt(DL - Dam_L - Gate_width, Dam_H + BW);			/**< Left top. */
Vec2d ConstrainP_rt(DL - Dam_L, Dam_H + BW);						/**< Right top. */
Vec2d ConstrainP_rb(DL - Dam_L, Base_bottom_position);				/**< Right bottom. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;						   /**< Reference density of fluid. */
Real gravity_g = 9.8e-3;				   /**< Value of gravity. */
Real U_f = 1.0;							   /**< Characteristic velocity. */
Real c_f = 20.0 * sqrt(140.0 * gravity_g); /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Material parameters of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 1.1;	 /**< Reference density of gate. */
Real poisson = 0.47; /**< Poisson ratio. */
Real Ae = 7.8e3;	 /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//----------------------------------------------------------------------
//	Fluid body with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_block_shape;
		water_block_shape.push_back(DamP_lb);
		water_block_shape.push_back(DamP_lt);
		water_block_shape.push_back(DamP_rt);
		water_block_shape.push_back(DamP_rb);
		water_block_shape.push_back(DamP_lb);
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Wall boundary body cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, -BW));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	create a gate shape
//----------------------------------------------------------------------
std::vector<Vecd> createGateShape()
{
	std::vector<Vecd> gate_shape;
	gate_shape.push_back(GateP_lb);
	gate_shape.push_back(GateP_lt);
	gate_shape.push_back(GateP_rt);
	gate_shape.push_back(GateP_rb);
	gate_shape.push_back(GateP_lb);

	return gate_shape;
}
//----------------------------------------------------------------------
//	Define the elastic gate body.
//----------------------------------------------------------------------
class Gate : public SolidBody
{
public:
	Gate(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 2.0))
	{
		/** Geometry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createGateShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
// Create the gate constrain shape
//----------------------------------------------------------------------
MultiPolygon createGateConstrainShape()
{
	//geometry
	std::vector<Vecd> gate_constrain_shape;
	gate_constrain_shape.push_back(ConstrainP_lb);
	gate_constrain_shape.push_back(ConstrainP_lt);
	gate_constrain_shape.push_back(ConstrainP_rt);
	gate_constrain_shape.push_back(ConstrainP_rb);
	gate_constrain_shape.push_back(ConstrainP_lb);

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(gate_constrain_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
//----------------------------------------------------------------------
//	Observer particle generator.
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(GateP_lb, 0.0));
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Set the starting time to zero. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** I/O environment. */
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	WaterBlock water_block(system, "WaterBody");
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f));

	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_boundary_particles(wall_boundary);

	Gate gate(system, "Gate");
	ElasticSolidParticles gate_particles(gate, makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson));
	/** offset particle position */
	gate_particles.offsetInitialParticlePosition(offset);

	ObserverBody gate_observer(system, "Observer", makeShared<SPHAdaptation>(1.15, 2.0));
	ObserverParticles observer_particles(gate_observer, makeShared<ObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexBodyRelation water_block_complex_relation(water_block, {&wall_boundary, &gate});
	BodyRelationInner gate_inner_relation(gate);
	BodyRelationContact gate_water_contact_relation(gate, {&water_block});
	BodyRelationContact gate_observer_contact_relation(gate_observer, {&gate});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Define the external force. */
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block, gravity);
	//----------------------------------------------------------------------
	//	Algorithms of fluid dynamics.
	//----------------------------------------------------------------------
	/** Evaluation of fluid density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex_relation);
	/** Compute time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Compute time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex_relation);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex_relation);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	/** Corrected configuration for solid dynamics. */
	solid_dynamics::CorrectConfiguration gate_corrected_configuration(gate_inner_relation);
	/** Compute the force exerted on elastic gate due to fluid pressure. */
	solid_dynamics::FluidPressureForceOnSolidRiemann fluid_pressure_force_on_gate(gate_water_contact_relation);
	//----------------------------------------------------------------------
	//	Algorithms of Elastic dynamics.
	//----------------------------------------------------------------------
	/** Compute time step size of elastic solid. */
	solid_dynamics::AcousticTimeStepSize gate_computing_time_step_size(gate);
	/** Stress relaxation stepping for the elastic gate. */
	solid_dynamics::StressRelaxationFirstHalf gate_stress_relaxation_first_half(gate_inner_relation);
	solid_dynamics::StressRelaxationSecondHalf gate_stress_relaxation_second_half(gate_inner_relation);
	/**Constrain a solid body part.  */
	MultiPolygonShape gate_constrain_shape(createGateConstrainShape());
	BodyRegionByParticle gate_constrain_part(gate, "GateConstrain", gate_constrain_shape);
	solid_dynamics::ConstrainSolidBodyRegion gate_constrain(gate, gate_constrain_part);
	/** Update the surface normal direaction of elastic gate. */
	solid_dynamics::UpdateElasticNormalDirection gate_update_normal(gate);
	/** Compute the average velocity of gate. */
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output body states for visualization in Tecplot. */
	BodyStatesRecordingToPlt write_real_body_states_to_plt(in_output, system.real_bodies_);
	/** Output body states for visualization in Paraview. */
	BodyStatesRecordingToVtp write_real_body_states_to_vtp(in_output, system.real_bodies_);
	/** Output the observed displacement of gate free end. */
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<indexVector, Vecd>>
		write_beam_tip_displacement("Position", in_output, gate_observer_contact_relation);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_particles.initializeNormalDirectionFromBodyShape();
	gate_particles.initializeNormalDirectionFromBodyShape();
	gate_corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int number_of_iterations = 0;
	int screen_output_interval = 100;
	Real End_Time = 400.0;			/**< End time. */
	Real D_Time = End_Time / 200.0; /**< time stamps for output. */
	Real dt = 0.0;					/**< Default acoustic time step sizes. */
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_real_body_states_to_vtp.writeToFile();
	write_beam_tip_displacement.writeToFile();
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
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			/** Update normal direction at elastic body surface. */
			gate_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** Fluid relaxation and force computation. */
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_gate.parallel_exec();
				density_relaxation.parallel_exec(dt);
				/** Solid dynamics time stepping. */
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					dt_s = gate_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s)
						dt_s = dt - dt_s_sum;
					gate_stress_relaxation_first_half.parallel_exec(dt_s);
					gate_constrain.parallel_exec();
					gate_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedList();
			gate.updateCellLinkedList();
			water_block_complex_relation.updateConfiguration();
			gate_water_contact_relation.updateConfiguration();
			/** Output the observed data. */
			write_beam_tip_displacement.writeToFile(number_of_iterations);
		}
		tick_count t2 = tick_count::now();
		write_real_body_states_to_vtp.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_beam_tip_displacement.newResultTest();

	return 0;
}
