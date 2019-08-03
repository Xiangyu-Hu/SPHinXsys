/**
 * @file 	elastic_gate.cpp
 * @brief 	2D elastic gate deformation due to dam break force.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid-structure-interaction (FSI) similation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
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
Real DL = 500.0; 									/**< Tank length. */
Real DH = 200.1; 									/**< Tank height. */
Real Dam_L = 100.0; 								/**< Dam width. */
Real Dam_H = 140.0; 								/**< Dam height. */
Real Rubber_width = 5.0;							/**< Width of the rubber gate. */
Real Base_bottom_position = 79.0;					/**< Position of gate base. (In Y direction) */
Real particle_spacing_ref = Rubber_width / 2.0; 	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4.0; 				/**< Extending width for BCs. */
/** The offset that the rubber gate shifted above the tank. */
Real dp_s = 0.5 * particle_spacing_ref;
Vec2d offset = Vec2d(0.0, Base_bottom_position - floor(Base_bottom_position / dp_s) * dp_s);
/**
 * @brief 	Define the corner point of dam geomerty.
 */
Vec2d DamP_lb(DL - Dam_L, 0.0); 		/**< Left bottom. */
Vec2d DamP_lt(DL - Dam_L, Dam_H); 		/**< Left top. */
Vec2d DamP_rt(DL, Dam_H); 				/**< Right top. */
Vec2d DamP_rb(DL, 0.0); 				/**< Right bottom. */
/**
 * @brief 	Define the corner point of gate base geomerty.
 */
Vec2d BaseP_lb(DL - Dam_L - Rubber_width, Base_bottom_position); 	/**< Left bottom. */
Vec2d BaseP_lt(DL - Dam_L - Rubber_width, DH); 						/**< Left top. */
Vec2d BaseP_rt(DL - Dam_L, DH); 									/**< Right top. */
Vec2d BaseP_rb(DL - Dam_L, Base_bottom_position); 					/**< Right bottom. */
/**
 * @brief 	Define the corner point of gate geomerty.
 */
Vec2d GateP_lb(DL - Dam_L - Rubber_width, 0.0); 					/**< Left bottom. */
Vec2d GateP_lt(DL - Dam_L - Rubber_width, Base_bottom_position); 	/**< Left top. */
Vec2d GateP_rt(DL - Dam_L, Base_bottom_position); 					/**< Right top. */
Vec2d GateP_rb(DL - Dam_L, 0.0); 									/**< Right bottom. */

/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;							/**< Reference density of fluid. */
Real gravity_g = 9.8e-3; 					/**< value of gravity. */
Real U_f = 1.0;								/**< Characteristic velocity. */
Real c_f = 20.0*sqrt(140.0*gravity_g); 		/**< Reference sound speed. */
Real mu_f = 0.0;							/**< Dynamics viscosity. */
Real k_f = 0.0;								/**< Thermal conduction rate. */

Real initial_pressure = 0.0;			/**< Initial pressure field. */
Vec2d intial_velocity(0.0, 0.0);		/**< Initial velocity field. */
/**
 * @brief Material properties of the elastic gate.
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
	WaterBlock(SPHSystem &system, string body_name,
		WeaklyCompressibleFluid &fluid_material, 
		int refinement_level, ParticlesGeneratorOps op)
		: FluidBody(system, body_name, fluid_material, refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> water_block_shape;
		water_block_shape.push_back(DamP_lb);
		water_block_shape.push_back(DamP_lt);
		water_block_shape.push_back(DamP_rt);
		water_block_shape.push_back(DamP_rb);
		water_block_shape.push_back(DamP_lb);
		Geometry *water_block_geometry = new Geometry(water_block_shape);
		body_region_.add_geometry(water_block_geometry, RegionBooleanOps::add);

		/** Finish the region modeling. */
		body_region_.done_modeling();
	}
};
/**
 * @brief 	wall body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, *(new Solid("EmptyWallMaterial")), refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> outer_wall_shape;
		outer_wall_shape.push_back(Point(-BW, -BW));
		outer_wall_shape.push_back(Point(-BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, DH + BW));
		outer_wall_shape.push_back(Point(DL + BW, -BW));
		outer_wall_shape.push_back(Point(-BW, -BW));
		body_region_.add_geometry(new Geometry(outer_wall_shape), RegionBooleanOps::add);

		std::vector<Point> inner_wall_shape;
		inner_wall_shape.push_back(Point(0.0, 0.0));
		inner_wall_shape.push_back(Point(0.0, DH));
		inner_wall_shape.push_back(Point(DL, DH));
		inner_wall_shape.push_back(Point(DL, 0.0));
		inner_wall_shape.push_back(Point(0.0, 0.0));
		body_region_.add_geometry(new Geometry(inner_wall_shape), RegionBooleanOps::sub);
		/** Finish the region modeling. */
		body_region_.done_modeling();

	}
};
/**
 * @brief  Gate base body definition.
 */
class GateBase : public SolidBody
{
public:
	GateBase(SPHSystem &system, string body_name, ElasticSolid &solid_material,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, solid_material, refinement_level, op)
	{
		/** Geometry definition. */
		std::vector<Point> gate_base_shape;
		gate_base_shape.push_back(BaseP_lb);
		gate_base_shape.push_back(BaseP_lt);
		gate_base_shape.push_back(BaseP_rt);
		gate_base_shape.push_back(BaseP_rb);
		gate_base_shape.push_back(BaseP_lb);
		body_region_.add_geometry(new Geometry(gate_base_shape), RegionBooleanOps::add);
		/** Finish the region modeling. */
		body_region_.done_modeling();
	}
};
/**
 * @brief  Define the elastic gate body.
 */
class Gate : public SolidBody
{
public:
	Gate(SPHSystem &system, string body_name, ElasticSolid &solid_material,
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, solid_material, refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> gate_shape;
		gate_shape.push_back(GateP_lb);
		gate_shape.push_back(GateP_lt);
		gate_shape.push_back(GateP_rt);
		gate_shape.push_back(GateP_rb);
		gate_shape.push_back(GateP_lb);
		body_region_.add_geometry(new Geometry(gate_shape), RegionBooleanOps::add);
		/** Finish the region modeling. */
		body_region_.done_modeling();
	}
};
/**
 * @brief Define the observer body.
 */
class Observer : public ObserverLagrangianBody
{
public:
	Observer(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: ObserverLagrangianBody(system, body_name, refinement_level, op)
	{
		/** Add observation point. */
		body_input_points_volumes_.push_back(make_pair(GateP_lb, 0.0));

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
	/** Define the external force. */
	Gravity gravity(Vecd(0.0, -gravity_g));
	/**
	 * @brief Material property, partilces and body creation of fluid.
	 */
	WeaklyCompressibleFluid 			fluid("Water", rho0_f, c_f, mu_f, k_f);
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", 
		fluid, 0, ParticlesGeneratorOps::lattice);
	FluidParticles 	fluid_particles(water_block);
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	WallBoundary *wall_boundary 
		= new WallBoundary(system, "Wall", 0, ParticlesGeneratorOps::lattice);
	SolidParticles 					solid_particles(wall_boundary);
	/**
	 * @brief 	Material property, particle and body creation of gate base.
	 */
	ElasticSolid 			solid_material("ElasticSolid", rho0_s, Youngs_modulus, poisson);
	GateBase *gate_base = new GateBase(system, "GateBase", 
		solid_material, 1, ParticlesGeneratorOps::lattice);
	ElasticSolidParticles 	gate_base_particles(gate_base);
	/**
	 * @brief 	Material property, particle and body creation of elastic gate.
	 */
	Gate *gate = new Gate(system, "Gate", 
		solid_material, 1, ParticlesGeneratorOps::lattice);
	ElasticSolidParticles 	gate_particles(gate);
	/**
	 * @brief 	Particle and body creation of gate observer.
	 */
	Observer *gate_observer 
		= new Observer(system, "Observer", 1, ParticlesGeneratorOps::direct);
	ObserverParticles 			observer_particles(gate_observer);
	/**
	 * @brief 	Body contact map.
	 * @details The contact map gives the data conntections between the bodies.
	 * 			Basically the the rang of bidies to build neighbor particle lists.
	 */
	SPHBodyTopology body_topology = { { water_block, { wall_boundary, gate_base, gate } },
									  { wall_boundary, { } },{ gate_base, { gate } },
									  { gate, { gate_base, water_block} }, { gate_observer,{ gate } } };
	system.SetBodyTopology(&body_topology);
	/**
	 * @brief 	Simulation set up.
	 */
	system.SetupSPHSimulation();
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /**
	  * @brief 	Methods used only once.
	  */
	  /** initial condition for fluid body */
	fluid_dynamics::WeaklyCompressibleFluidInitialCondition set_all_fluid_particles_at_rest(water_block);
	/** obtain the initial number density of fluid. */
	fluid_dynamics::InitialNumberDensity 		fluid_initial_number_density(water_block,
		{ wall_boundary,  gate_base, gate });

	/** initial condition for the solid body */
	solid_dynamics::SolidDynamicsInitialCondition set_all_wall_particles_at_rest(wall_boundary);
	/** initial condition for the elastic solid bodies */
	solid_dynamics::ElasticSolidDynamicsInitialCondition set_all_gate_base_particles_at_rest(gate_base);
	solid_dynamics::ElasticSolidDynamicsInitialCondition set_all_gate_particles_at_rest(gate);
	/** Initialize normal direction of the wall boundary. */
	solid_dynamics::NormalDirectionSummation 	get_wall_normal(wall_boundary, {});
	/** Initialize normal direction of the gate base. */
	solid_dynamics::NormalDirectionSummation 	get_gate_base_normal(gate_base, { gate });
	/** Initialize normal direction of the elastic gate. */
	solid_dynamics::NormalDirectionSummation 	get_gate_normal(gate, { gate_base });
	/** Corrected strong configuration. */
	solid_dynamics::CorrectConfiguration 		gate_base_corrected_configuration_in_strong_form(gate_base, { gate });
	solid_dynamics::CorrectConfiguration 		gate_corrected_configuration_in_strong_form(gate, { gate_base });


	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	InitializeOtherAccelerations 	initialize_fluid_acceleration(water_block, &gravity);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of fluid density by summation approach. */
	fluid_dynamics::DensityBySummationFreeSurface		update_fluid_desnity(water_block, { wall_boundary,  gate_base, gate });
	/** Compute time step size without considering sound wave speed. */
	fluid_dynamics::GetAdvectionTimeStepSize			get_fluid_adevction_time_step_size(water_block, U_f);
	/** Compute time step size with considering sound wave speed. */
	fluid_dynamics::GetAcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationVerletFreeSurface pressure_relaxation(water_block,
		{ wall_boundary,  gate_base, gate }, &gravity);
	/**
	 * @brief Algorithms of FSI.
	 */
	 /** Compute the force exerted on elastic gate due to fluid pressure. */
	solid_dynamics::FluidPressureForceOnSolid 	fluid_pressure_force_on_gate(gate, { water_block }, &gravity);
	/**
	 * @brief Algorithms of Elastic dynamics.
	 */
	 /** Compute time step size of elastic solid. */
	solid_dynamics::GetAcousticTimeStepSize 	gate_computing_time_step_size(gate);
	/** Stress relaxation stepping for the elastic gate. */
	solid_dynamics::StressRelaxation 			gate_stress_relaxation(gate, { gate_base });
	/** Stress update for contrained wall body(gate base). */
	solid_dynamics::StressInConstrinedElasticBodyFirstHalf 	gate_base_stress_update_first_half(gate_base);
	solid_dynamics::StressInConstrinedElasticBodySecondHalf gate_base_stress_update_second_half(gate_base, { gate });
	/** Update the norm of elastic gate. */
	solid_dynamics::UpdateElasticNormalDirection 	gate_update_normal(gate);
	/** Compute the average velocity of gate. */
	solid_dynamics::InitializeDisplacement 			gate_initialize_displacement(gate);
	solid_dynamics::UpdateAverageVelocity 			gate_average_velocity(gate);
	/**
	 * @brief 	Methods used for updating data structure.
	 */
	 /** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList 				update_water_block_cell_linked_list(water_block);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsConfiguration 				update_water_block_configuration(water_block);
	/** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList 				update_gate_cell_linked_list(gate);
	/** Update the contact configuration for a given contact map. */
	ParticleDynamicsInteractionConfiguration 	update_gate_interaction_configuration(gate, { water_block });
	/**
	 * @brief Output.
	 */
	In_Output in_output(system);
	/** Output body states for visulaization. */
	WriteBodyStatesToPlt 				write_real_body_states_to_plt(in_output, system.real_bodies_);
	/** Output body states for visulaization. */
	WriteBodyStatesToVtu 				write_real_body_states_to_vtu(in_output, system.real_bodies_);
	/** Output the observed displacement of gate free end. */
	WriteObservedElasticDisplacement 	write_beam_tip_displacement(in_output, gate_observer, { gate });
	/**
	 * @brief The time stepping starts here.
	 */
	 /** Set the starting time to zero. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	set_all_fluid_particles_at_rest.exec();
	set_all_wall_particles_at_rest.exec();
	set_all_gate_base_particles_at_rest.exec();
	set_all_gate_particles_at_rest.exec();

	get_wall_normal.parallel_exec();
	get_gate_base_normal.parallel_exec();
	get_gate_normal.parallel_exec();
	fluid_initial_number_density.parallel_exec();
	gate_corrected_configuration_in_strong_form.parallel_exec();
	gate_base_corrected_configuration_in_strong_form.parallel_exec();
	/**
	 * @brief Initial output.
	 */
	write_real_body_states_to_plt.WriteToFile(GlobalStaticVariables::physical_time_);
	write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real End_Time = 400.0;			/**< End time. */
	Real D_Time = End_Time / 200.0;	/**< time stamps for output. */
	Real Dt = 0.0;					/**< Default advection time step sizes. */
	Real dt = 0.0; 					/**< Default accoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	/**
	 * @brief Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integeral_time < D_Time)
		{
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			/** Acceleration due to viscous force and gravity. */
			initialize_fluid_acceleration.parallel_exec();
			/** Update normal direction on elastic body. */
			gate_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}
				/** Fluid relaxation and force computaton. */
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_gate.parallel_exec();
				/** Solid dynamics time stepping. */
				Real dt_s_sum = 0.0;
				gate_initialize_displacement.parallel_exec();
				//if(GlobalStaticVariables::physical_time_ >= 100.0){
				while (dt_s_sum < dt)
				{
					Real dt_s = gate_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;
					if (ite % 100 == 0) {
						cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_ << "	dt_s: "
							<< dt_s << "\n";
					}
					gate_base_stress_update_first_half.parallel_exec(dt_s);
					gate_stress_relaxation.parallel_exec(dt_s);
					gate_base_stress_update_second_half.exec(dt_s);
					dt_s_sum += dt_s;
				}
				//}
				gate_average_velocity.parallel_exec(dt);

				ite++;
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			/** Update cell linked list and configuration. */
			update_water_block_cell_linked_list.parallel_exec();
			update_water_block_configuration.parallel_exec();
			update_gate_cell_linked_list.parallel_exec();
			update_gate_interaction_configuration.parallel_exec();
			/** Output the observed data. */
			write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);
		}
		tick_count t2 = tick_count::now();
		write_real_body_states_to_vtu.WriteToFile(GlobalStaticVariables::physical_time_  * 0.001);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
