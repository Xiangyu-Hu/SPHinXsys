/**
 * @file 	elastic_gate.cpp
 * @brief 	2D elastic gate deformation due to dam break force.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid-structure-interaction (FSI) simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
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
Real DL = 500.0; 									/**< Tank length. */
Real DH = 200.1; 									/**< Tank height. */
Real Dam_L = 100.0; 								/**< Water block width. */
Real Dam_H = 140.0; 								/**< Water block height. */
Real Gate_width = 5.0;								/**< Width of the gate. */
Real Base_bottom_position = 79.0;					/**< Position of gate base. (In Y direction) */
Real resolution_ref = Gate_width / 2.0; 			/**< Initial reference particle spacing. */
Real BW = resolution_ref * 4.0; 					/**< Extending width for BCs. */
/** The offset that the rubber gate shifted above the tank. */
Real dp_s = 0.5 * resolution_ref;
Vec2d offset = Vec2d(0.0, Base_bottom_position - floor(Base_bottom_position / dp_s) * dp_s);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/**
 * @brief 	Define the corner points of the water block geomtry.
 */
Vec2d DamP_lb(DL - Dam_L, 0.0); 		/**< Left bottom. */
Vec2d DamP_lt(DL - Dam_L, Dam_H); 		/**< Left top. */
Vec2d DamP_rt(DL, Dam_H); 				/**< Right top. */
Vec2d DamP_rb(DL, 0.0); 				/**< Right bottom. */
/**
 * @brief 	Define the corner points of the gate geomtry.
 */
Vec2d GateP_lb(DL - Dam_L - Gate_width, 0.0); 				/**< Left bottom. */
Vec2d GateP_lt(DL - Dam_L - Gate_width, Dam_H + BW); 		/**< Left top. */
Vec2d GateP_rt(DL - Dam_L, Dam_H + BW); 					/**< Right top. */
Vec2d GateP_rb(DL - Dam_L, 0.0); 							/**< Right bottom. */
/**
 * @brief 	Define the corner points of the gate constrain.
 */
Vec2d ConstrainP_lb(DL - Dam_L - Gate_width, Base_bottom_position); 	/**< Left bottom. */
Vec2d ConstrainP_lt(DL - Dam_L - Gate_width, Dam_H + BW); 				/**< Left top. */
Vec2d ConstrainP_rt(DL - Dam_L, Dam_H + BW); 							/**< Right top. */
Vec2d ConstrainP_rb(DL - Dam_L, Base_bottom_position); 					/**< Right bottom. */

/**
 * @brief Material parameters of the fluid.
 */
Real rho0_f = 1.0;							/**< Reference density of fluid. */
Real gravity_g = 9.8e-3; 					/**< Value of gravity. */
Real U_f = 1.0;								/**< Characteristic velocity. */
Real c_f = 20.0*sqrt(140.0*gravity_g); 		/**< Reference sound speed. */
/**
 * @brief Material parameters of the elastic gate.
 */
Real rho0_s = 1.1; 						/**< Reference density of gate. */
Real poisson = 0.47; 					/**< Poisson ratio. */
Real Ae = 7.8e3; 						/**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
/**
 * @brief 	Fluid body definition.
 */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, std::string body_name)
		: FluidBody(system, body_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> water_block_shape;
		water_block_shape.push_back(DamP_lb);
		water_block_shape.push_back(DamP_lt);
		water_block_shape.push_back(DamP_rt);
		water_block_shape.push_back(DamP_rb);
		water_block_shape.push_back(DamP_lb);
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
/**
 * @brief 	Case-dependent material properties definition.
 */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial()	: WeaklyCompressibleFluid()
	{
		rho0_ = rho0_f;
		c0_ = c_f;

		assignDerivedMaterialParameters();
	}
};
/**
 * @brief 	wall body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, std::string body_name) : SolidBody(system, body_name)
	{
		/** Geomtry definition. */
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

		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
/**
* @brief create a gate shape
*/
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
/**
 * @brief  Define the elastic gate body.
 */
class Gate : public SolidBody
{
public:
	Gate(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(1.15, 1))
	{
		/** Geomtry definition. */
		std::vector<Vecd> gate_shape = createGateShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(gate_shape, ShapeBooleanOps::add);
	}
};
/**
* @brief create the gate constrain shape
*/
std::vector<Vecd> createGateConstrainShape()
{
	//geometry
	std::vector<Vecd> gate_constrain_shape;
	gate_constrain_shape.push_back(ConstrainP_lb);
	gate_constrain_shape.push_back(ConstrainP_lt);
	gate_constrain_shape.push_back(ConstrainP_rt);
	gate_constrain_shape.push_back(ConstrainP_rb);
	gate_constrain_shape.push_back(ConstrainP_lb);

	return gate_constrain_shape;
}
/**
* @brief define the beam base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class GateConstrain : public BodyPartByParticle
{
public:
	GateConstrain(SolidBody* solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		/* Geometry definition */
		std::vector<Vecd> gate_constrain_shape = createGateConstrainShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(gate_constrain_shape, ShapeBooleanOps::add);

		//tag the constrained particle
		tagBodyPart();
	}
};
/**
 * @brief Define gate material.
 */
class GateMaterial : public LinearElasticSolid
{
public:
	GateMaterial() : LinearElasticSolid()
	{
		rho0_ = rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};
/**
 * @brief Define the observer body.
 */
class Observer : public FictitiousBody
{
public:
	Observer(SPHSystem &system, std::string body_name) :
		FictitiousBody(system, body_name, new ParticleAdaptation(1.15, 1))
	{
		/** Add observation point. */
		body_input_points_volumes_.push_back(std::make_pair(GateP_lb, 0.0));

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
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Define the external force. */
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** Set the starting time to zero. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	WaterBlock *water_block = new WaterBlock(system, "WaterBody");
	WaterMaterial 	*water_material = new WaterMaterial();
	FluidParticles 	fluid_particles(water_block, water_material);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall");
	SolidParticles 		wall_boundary_particles(wall_boundary);
	/**
	 * @brief 	Material property, particle and body creation of gate.
	 */
	GateMaterial* gate_material = new GateMaterial();
	Gate *gate = new Gate(system, "Gate");
	ElasticSolidParticles 	gate_particles(gate, gate_material);
	/** offset particle position */
	gate_particles.offsetInitialParticlePosition(offset);
	/**
	 * @brief 	Particle and body creation of gate observer.
	 */
	Observer *gate_observer = new Observer(system, "Observer");
	BaseParticles 			observer_particles(gate_observer);

	/** topology */
	ComplexBodyRelation* water_block_complex_relation = new ComplexBodyRelation(water_block, { wall_boundary, gate});
	BodyRelationInner*	gate_inner_relation = new BodyRelationInner(gate);
	BodyRelationContact* gate_water_contact_relation = new BodyRelationContact(gate, { water_block });
	BodyRelationContact* gate_observer_contact_relation = new BodyRelationContact(gate_observer, { gate });

	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /**
	  * @brief 	Methods used only once.
	  */
	/** Corrected strong configuration. */
	solid_dynamics::CorrectConfiguration 		gate_corrected_configuration_in_strong_form(gate_inner_relation);

	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	TimeStepInitialization 	initialize_a_fluid_step(water_block, &gravity);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of fluid density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex	update_density_by_summation(water_block_complex_relation);
	/** Compute time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize			get_fluid_advection_time_step_size(water_block, U_f);
	/** Compute time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex_relation);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex_relation);
	/**
	 * @brief Algorithms of FSI.
	 */
	 /** Compute the force exerted on elastic gate due to fluid pressure. */
	solid_dynamics::FluidPressureForceOnSolidRiemann 	fluid_pressure_force_on_gate(gate_water_contact_relation);
	/**
	 * @brief Algorithms of Elastic dynamics.
	 */
	 /** Compute time step size of elastic solid. */
	solid_dynamics::AcousticTimeStepSize 	gate_computing_time_step_size(gate);
	/** Stress relaxation stepping for the elastic gate. */
	solid_dynamics::StressRelaxationFirstHalf 			gate_stress_relaxation_first_half(gate_inner_relation);
	solid_dynamics::StressRelaxationSecondHalf 			gate_stress_relaxation_second_half(gate_inner_relation);
	/**Constrain a solid body part.  */
	solid_dynamics::ConstrainSolidBodyRegion
		gate_constrain(gate, new GateConstrain(gate, "GateConstrain"));
	/** Update the norm of elastic gate. */
	solid_dynamics::UpdateElasticNormalDirection 	gate_update_normal(gate);
	/** Compute the average velocity of gate. */
	solid_dynamics::AverageVelocityAndAcceleration 		average_velocity_and_acceleration(gate);

	/**
	 * @brief Output.
	 */
	In_Output in_output(system);
	/** Output body states for visualization. */
	BodyStatesRecordingToPlt 				write_real_body_states_to_plt(in_output, system.real_bodies_);
	/** Output body states for visualization. */
	BodyStatesRecordingToVtu 				write_real_body_states_to_vtu(in_output, system.real_bodies_);
	/** Output the observed displacement of gate free end. */
	ObservedQuantityRecording<indexVector, Vecd>
		write_beam_tip_displacement("Position", in_output, gate_observer_contact_relation);
	/**
	 * @brief The time stepping starts here.
	 */
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_particles.initializeNormalDirectionFromGeometry();
	gate_particles.initializeNormalDirectionFromGeometry();
	gate_corrected_configuration_in_strong_form.parallel_exec();

	write_real_body_states_to_vtu.writeToFile(0);
	write_beam_tip_displacement.writeToFile(0);

	int number_of_iterations = 0;
	int screen_output_interval = 100;
	Real End_Time = 400.0;			/**< End time. */
	Real D_Time = End_Time / 200.0;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes. */
	Real dt = 0.0; 					/**< Default acoustic time step sizes. */
	Real dt_s = 0.0;				/**< Default acoustic time step sizes for solid. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	/**
	 * @brief Main loop starts here.
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
			update_density_by_summation.parallel_exec();
			/** Update normal direction on elastic body. */
			gate_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** Fluid relaxation and force computaton. */
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_gate.parallel_exec();
				density_relaxation.parallel_exec(dt);
				/** Solid dynamics time stepping. */
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					dt_s = gate_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;
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
			water_block->updateCellLinkedList();
			gate->updateCellLinkedList();
			water_block_complex_relation->updateConfiguration();
			gate_water_contact_relation->updateConfiguration();
			/** Output the observed data. */
			write_beam_tip_displacement.writeToFile(number_of_iterations);
		}
		tick_count t2 = tick_count::now();
		write_real_body_states_to_vtu.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
