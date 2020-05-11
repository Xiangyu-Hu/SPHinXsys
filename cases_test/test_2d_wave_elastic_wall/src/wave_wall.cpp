/**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
Real DL = 18.42; //tank length
Real DH = 1.0; //tank height
Real DL_Extra = 1.0; // for wave maker

Real Water_H = 0.69; //dam height

Real Flap_x     = 7.92;
Real Flap_H = 0.48;
Real Flap_width = 0.12;

Real Base_bottom_position = 0.15;
Real Base_constrain_height = 0.1;
Real particle_spacing_ref = Flap_width / 4.0; //particle spacing
Real BW = particle_spacing_ref * 4.0; //boundary width

//the offset that the rubber gate shifted above the tank
Vec2d offset = Vec2d(0.0, Base_bottom_position 
	- floor(Base_bottom_position / particle_spacing_ref) * particle_spacing_ref);

//define Dam domain
Vec2d Water_lb(0.0, 0.0); //left bottom
Vec2d Water_lt(0.0, Water_H); //left top
Vec2d Water_rt(DL, Water_H); //right top
Vec2d Water_rb(DL, 0.35); //right bottom
Vec2d Water_slope_1(DL - 6.2, 0.35);
Vec2d Water_slope_2(DL - 6.2 - 3.7, 0.35 - 0.2);
Vec2d Water_slope_3(DL - 6.2 - 3.7 - 2.4, 0.35 - 0.2);
Vec2d Water_slope_4(DL - 6.2 - 3.7 - 2.4 - 1.3, 0.0);

//define constrin region of gate
Vec2d Base_lb(Flap_x - 0.5 * Flap_width, Base_bottom_position); //left bottom
Vec2d Base_lt(Flap_x - 0.5 * Flap_width, Base_bottom_position + Base_constrain_height); //left top
Vec2d Base_rt(Flap_x + 0.5 * Flap_width, Base_bottom_position + Base_constrain_height); //right top
Vec2d Base_rb(Flap_x + 0.5 * Flap_width, Base_bottom_position ); //right bottom

//define gate
Vec2d Flap_lb(Flap_x - 0.5 * Flap_width, Base_bottom_position); //left bottom
Vec2d Flap_lt(Flap_x - 0.5 * Flap_width, Base_bottom_position + Base_constrain_height + Flap_H ); //left top
Vec2d Flap_rt(Flap_x + 0.5 * Flap_width, Base_bottom_position + Base_constrain_height + Flap_H ); //right top
Vec2d Flap_rb(Flap_x + 0.5 * Flap_width, Base_bottom_position ); //right bottom

//gravity value
Real gravity_g = 9.8;

//for material properties of the fluid
Real rho0_f = 1000.0;
Real U_f = 2.0*sqrt(0.79 * gravity_g);
Real c_f = 10.0*U_f;

//for material properties of the solid
Real rho0_s = 1100.0; //reference density
Real poisson = 0.47; //Poisson ratio
Real Youngs_modulus = 7.8e6; //normalized Youngs Modulus
/**
* @brief define geometry and initial conditions of SPH bodies
*/
/**
* @brief create a water block shape
*/
std::vector<Point> CreatWaterBlockShape()
{
	//geometry
	std::vector<Point> pnts;
	pnts.push_back(Water_lb);
	pnts.push_back(Water_lt);
	pnts.push_back(Water_rt);
	pnts.push_back(Water_rb);
	pnts.push_back(Water_slope_1);
	pnts.push_back(Water_slope_2);
	pnts.push_back(Water_slope_3);
	pnts.push_back(Water_slope_4);
	pnts.push_back(Water_lb);

	return pnts;
}
/**
* @brief create gate base shape
*/
std::vector<Point> CreatGateConstrainShape()
{
	//geometry
	std::vector<Point> pnts2;
	pnts2.push_back(Base_lb);
	pnts2.push_back(Base_lt);
	pnts2.push_back(Base_rt);
	pnts2.push_back(Base_rb);
	pnts2.push_back(Base_lb);

	return pnts2;
}
/**
* @brief create gate base shape
*/
std::vector<Point> CreatFlapShape()
{
	//geometry
	std::vector<Point> pnts3;
	pnts3.push_back(Flap_lb);
	pnts3.push_back(Flap_lt);
	pnts3.push_back(Flap_rt);
	pnts3.push_back(Flap_rb);
	pnts3.push_back(Flap_lb);

	return pnts3;
}
/**
* @brief create outer wall shape
*/
std::vector<Point> CreatOuterWallShape()
{
	std::vector<Point> pnts1;
	pnts1.push_back(Point(-DL_Extra - BW, -BW));
	pnts1.push_back(Point(-DL_Extra - BW, DH + BW));
	pnts1.push_back(Point(DL + BW, DH + BW));
	pnts1.push_back(Point(DL + BW, 0.35 - BW));
	pnts1.push_back(Water_slope_1 + (0.0, -BW));
	pnts1.push_back(Water_slope_2 + (0.0, -BW));
	pnts1.push_back(Water_slope_3 + (0.0, -BW));
	pnts1.push_back(Water_slope_4 + (0.0, -BW));
	pnts1.push_back(Point(-DL_Extra - BW, -BW));

	return pnts1;
}
/**
* @brief create inner wall shape 01
*/
std::vector<Point> CreatInnerWallShape01()
{
	std::vector<Point> pnts2;
	pnts2.push_back(Water_lb);
	pnts2.push_back(Point(0.0, DH + BW));
	pnts2.push_back(Point(DL, DH + BW));
	pnts2.push_back(Water_rb);
	pnts2.push_back(Water_slope_1);
	pnts2.push_back(Water_slope_2);
	pnts2.push_back(Water_slope_3);
	pnts2.push_back(Water_slope_4);
	pnts2.push_back(Water_lb);

	return pnts2;
}
/**
* @brief create inner wall shape 02
*/
std::vector<Point> CreatInnerWallShape02()
{
	std::vector<Point> pnts3;
	pnts3.push_back(Point(-DL_Extra, 0.0));
	pnts3.push_back(Point(-DL_Extra, DH + BW));
	pnts3.push_back(Point(-BW, DH + BW));
	pnts3.push_back(Point(-BW, 0.0));
	pnts3.push_back(Point(-DL_Extra, 0.0));

	return pnts3;
}
/**
* @brief create wave maker shape
*/
std::vector<Point> CreatWaveMakerShape()
{
	std::vector<Point> wave_make_shape;
	wave_make_shape.push_back(Point(-BW, 0.0));
	wave_make_shape.push_back(Point(-BW, DH + BW));
	wave_make_shape.push_back(Point(0.0, DH + BW));
	wave_make_shape.push_back(Point(0.0, 0.0));
	wave_make_shape.push_back(Point(-BW, 0.0));

	return wave_make_shape;
}
//define the fluid body
class WaterBlock : public FluidBody
{
	public:
		WaterBlock(SPHSystem &system, string body_name,
			int refinement_level, ParticlesGeneratorOps op)
			: FluidBody(system, body_name, refinement_level, op)
		{
			std::vector<Point> water_block_shape = CreatWaterBlockShape();
			std::vector<Point> flap_shape = CreatFlapShape();
			body_region_.add_polygon(water_block_shape, RegionBooleanOps::add);
			body_region_.add_polygon(flap_shape, RegionBooleanOps::sub);
			//finish the region modeling
			body_region_.done_modeling();
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
//define the static solid wall boudary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> outer_wall_shape   = CreatOuterWallShape();
		std::vector<Point> inner_wall_shape_1 = CreatInnerWallShape01();
		std::vector<Point> inner_wall_shape_2 = CreatInnerWallShape02();
		body_region_.add_polygon(outer_wall_shape, RegionBooleanOps::add);
		body_region_.add_polygon(inner_wall_shape_1, RegionBooleanOps::sub);
		body_region_.add_polygon(inner_wall_shape_2, RegionBooleanOps::sub);
		//finish the region modeling
		body_region_.done_modeling();
	}
};

//define the elastic gate
class Gate : public SolidBody
{
public:
	Gate(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> flap_shape = CreatFlapShape();
		body_region_.add_polygon(flap_shape, RegionBooleanOps::add);
		//finish the region modeling
		body_region_.done_modeling();
	}
};

/**
* @brief define the beam base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class GateConstrain : public BodyPartByParticle
{
public:
	GateConstrain(SolidBody* solid_body, string constrianed_region_name)
		: BodyPartByParticle(solid_body, constrianed_region_name)
	{
		/* Geometry defination */
		std::vector<Point> beam_base_shape = CreatGateConstrainShape();
		Geometry* beam_base_gemetry = new Geometry(beam_base_shape);
		body_part_region_.add_geometry(beam_base_gemetry, RegionBooleanOps::add);
		/** Finish the region modeling. */
		body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
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
		rho_0_ = rho0_s;
		E_0_ = Youngs_modulus;
		nu_ = poisson;

		assignDerivedMaterialParameters();
	}
};
/**
* @brief constrain the beam base
*/
class WaveMaker : public BodyPartByParticle
{
public:
	WaveMaker(SolidBody *solid_body, string constrianed_region_name)
		: BodyPartByParticle(solid_body, constrianed_region_name)
	{
		//geometry
		std::vector<Point> wave_maker_shape = CreatWaveMakerShape();
		body_part_region_.add_polygon(wave_maker_shape, RegionBooleanOps::add);
		//finish the region modeling
		body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
	}
};
/**
* @brief making the wave
*/
class WaveMaking : public solid_dynamics::ConstrainSolidBodyRegion
{
	Real wave_freq_;
	Real wave_stroke_;
	Real time_;

	virtual Vec2d GetDisplacement(Vec2d &pos) override
	{
		Vec2d displacement(0);
		displacement[0] = 0.5 * wave_stroke_ 
			* sin(wave_freq_ * time_);
		return displacement;
	}

	virtual Vec2d GetVelocity(Vec2d &pos) override
	{
		Vec2d velcoity(0);
		velcoity[0] = 0.5 * wave_stroke_ * wave_freq_ 
			* cos(wave_freq_ * time_);
		return velcoity;
	}

	virtual Vec2d GetAcceleration(Vec2d &pos) override
	{
		Vec2d acceleration(0);
		acceleration[0] = -0.5 * wave_stroke_ * wave_freq_ * wave_freq_ 
			* sin(wave_freq_ * time_);
		return acceleration;
	}

	virtual void PrepareConstraint() override 
	{
		time_ = GlobalStaticVariables::physical_time_;
	}

public:
	WaveMaking(SolidBody *solid_body, BodyPartByParticle*constrained_region)
		: ConstrainSolidBodyRegion(solid_body, constrained_region), time_(0.0)
	{
		wave_freq_ = 3.14159;
		wave_stroke_ = 0.21515;
	}
};

//the main program
int main()
{
	//build up context -- a SPHSystem
	SPHSystem system(Vec2d( - DL - BW , -BW), 
		Vec2d(DL + BW, DH + BW), particle_spacing_ref);

	//define external force
	Gravity gravity(Vecd(0.0, -gravity_g));

	//the water block
	WaterBlock *water_block 
		= new WaterBlock(system, "WaterBody", 0, ParticlesGeneratorOps::lattice);
	//fluid material properties
	WaterMaterial *water_material = new WaterMaterial();
	//creat fluid particles
	FluidParticles fluid_particles(water_block, water_material);

	//the wall boundary
	WallBoundary *wall_boundary 
		= new WallBoundary(system, "Wall", 0, ParticlesGeneratorOps::lattice);
	//creat solid particles
	SolidParticles solid_particles(wall_boundary);

	//the elastic gate
	//elastic soild material properties
	GateMaterial* gate_material = new GateMaterial();
	Gate *gate =
		new Gate(system, "Gate", 0, ParticlesGeneratorOps::lattice);
	//creat particles for the elastic gate
	ElasticSolidParticles gate_particles(gate, gate_material);

	//set body contact map
	//the contact map gives the data conntections between the bodies
	//basically the the rang of bidies to build neighbor particle lists
	SPHBodyTopology body_topology 
		= { { water_block, { wall_boundary, gate } },
		{ wall_boundary, { gate} },
	    { gate, {wall_boundary, water_block} } };
	system.SetBodyTopology(&body_topology);

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	//methods only used only once
	//-------------------------------------------------------------------

	//initialize normal direction of the wall boundary
	solid_dynamics::NormalDirectionSummation get_wall_normal(wall_boundary, { gate });
	//initialize normal direction of the elastic gate
	solid_dynamics::NormalDirectionSummation get_gate_normal(gate, { wall_boundary });
	//corrected strong configuration	
	solid_dynamics::CorrectConfiguration
		gate_corrected_configuration_in_strong_form(gate);


	//--------------------------------------------------------------------------
	//methods used for time stepping
	//--------------------------------------------------------------------------
	/** Add particle acceleration due to gravity force. */
	//fluid dynamics
	InitializeATimeStep 	initialize_a_fluid_step(water_block, &gravity);
	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummationFreeSurface
		update_fluid_desnity(water_block, { wall_boundary,  gate });
	//time step size without considering sound wave speed
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationFirstHalfRiemann
		pressure_relaxation_first_half(water_block, { wall_boundary,  gate });
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(water_block, { wall_boundary, gate });

	//FSI
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_gate(gate, { water_block });

	//solid dynmaics
	//time step size caclutation
	solid_dynamics::GetAcousticTimeStepSize gate_computing_time_step_size(gate);
	//stress relaxation for the beam
	solid_dynamics::StressRelaxationFirstHalf gate_stress_relaxation_first_half(gate);
	solid_dynamics::StressRelaxationSecondHalf gate_stress_relaxation_second_half(gate);

	//average velocity
	solid_dynamics::InitializeDisplacement 		gate_initialize_displacement(gate);
	solid_dynamics::UpdateAverageVelocity 		gate_average_velocity(gate);
	solid_dynamics::UpdateElasticNormalDirection 	gate_update_normal(gate);

	/** Constrain the gate base.  */
	solid_dynamics::ConstrainSolidBodyRegion
		gate_constrain(gate, new GateConstrain(gate, "GateConstrain"));
	//constrain region of the part of wall boundary
	WaveMaking
		wave_making(wall_boundary, new WaveMaker(wall_boundary, "WaveMaker"));

	//-------------------------------------------------------------------
	//methods used for updating data structure
	//-------------------------------------------------------------------
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList update_water_block_cell_linked_list(water_block);
	//update the configuration of bodies when neccessary
	ParticleDynamicsConfiguration update_water_block_configuration(water_block);
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList update_wall_boundry_cell_linked_list(wall_boundary);
	//update the configuration of bodies when neccessary
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList update_gate_cell_linked_list(gate);
	//update the contact configuration for a given contact map
	ParticleDynamicsInteractionConfiguration
		update_gate_interaction_configuration(gate, { water_block });

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	WriteBodyStatesToVtu write_real_body_states(in_output, system.real_bodies_);

	//-------------------------------------------------------------------
	//from here the time stepping begines
	//-------------------------------------------------------------------
	//starting time zero
	GlobalStaticVariables::physical_time_ = 0.0;

	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	system.InitializeSystemCellLinkedLists();
	system.InitializeSystemConfigurations();
	get_wall_normal.parallel_exec();
	get_gate_normal.parallel_exec();
	gate_corrected_configuration_in_strong_form.parallel_exec();

	//initial output
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);

	int number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	Real End_Time = 10.0;
	//time step size for oupt file
	Real D_Time = End_Time/100.0;
	Real Dt = 0.0;//default advection time step sizes
	Real dt = 0.0; //default accoustic time step sizes
	Real dt_s = 0.0;	/**< Default acoustic time step sizes for solid. */

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < D_Time) {

			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			gate_update_normal.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				wave_making.parallel_exec(dt);
				//fluid dynamic, first half 
				pressure_relaxation_first_half.parallel_exec(dt);
				//FSI on pressure force
				fluid_pressure_force_on_gate.parallel_exec();
				//fluid dynamic, second half 
				pressure_relaxation_second_half.parallel_exec(dt);

				//solid dynamics
				Real dt_s_sum = 0.0;
				gate_initialize_displacement.parallel_exec();
				while (dt_s_sum < dt) 
				{
					dt_s = gate_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;
					gate_stress_relaxation_first_half.parallel_exec(dt_s);
					gate_constrain.parallel_exec();
					gate_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				gate_average_velocity.parallel_exec(dt);

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			
			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
			}
			number_of_iterations++;

			//water block confifuration
			update_water_block_cell_linked_list.parallel_exec();
			update_water_block_configuration.parallel_exec();

			//gate cell linked list and  interaction configuration
			update_gate_cell_linked_list.parallel_exec();
			update_gate_interaction_configuration.parallel_exec();

			update_wall_boundry_cell_linked_list.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_  * 0.001);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
	
	return 0;
}
