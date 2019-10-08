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

//define rigid base for installation of gate
Vec2d Base_lb(Flap_x - 0.5 * Flap_width, Base_bottom_position); //left bottom
Vec2d Base_lt(Flap_x - 0.5 * Flap_width, Base_bottom_position + 0.1); //left top
Vec2d Base_rt(Flap_x + 0.5 * Flap_width, Base_bottom_position + 0.1); //right top
Vec2d Base_rb(Flap_x + 0.5 * Flap_width, Base_bottom_position ); //right bottom

//define gate
Vec2d Flap_lb(Flap_x - 0.5 * Flap_width, Base_bottom_position + 0.1); //left bottom
Vec2d Flap_lt(Flap_x - 0.5 * Flap_width, Base_bottom_position + 0.1 + Flap_H ); //left top
Vec2d Flap_rt(Flap_x + 0.5 * Flap_width, Base_bottom_position + 0.1 + Flap_H ); //right top
Vec2d Flap_rb(Flap_x + 0.5 * Flap_width, Base_bottom_position + 0.1 ); //right bottom

//gravity value
Real gravity_g = 9.8;

//for material properties of the fluid
Real rho0_f = 1000.0;
Real U_f = 2.0*sqrt(0.79 * gravity_g);
Real c_f = 10.0*U_f;
Real mu_f = 0.0;
Real k_f = 0.0;

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
std::vector<Point> CreatGateBaseShape()
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
			std::vector<Point> gate_base_shape = CreatGateBaseShape();
			std::vector<Point> flap_shape = CreatFlapShape();
			body_region_.add_polygon(water_block_shape, RegionBooleanOps::add);
			body_region_.add_polygon(gate_base_shape, RegionBooleanOps::sub);
			body_region_.add_polygon(flap_shape, RegionBooleanOps::sub);
			//finish the region modeling
			body_region_.done_modeling();
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

//define base of the elastic gate
class GateBase : public SolidBody
{
public:
	GateBase(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> gate_base_shape = CreatGateBaseShape();
		body_region_.add_polygon(gate_base_shape, RegionBooleanOps::add);
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
* @brief constrain the beam base
*/
class WaveMaker : public SolidBodyPart
{
public:
	WaveMaker(SolidBody *solid_body, string constrianed_region_name)
		: SolidBodyPart(solid_body, constrianed_region_name)
	{
		//geometry
		std::vector<Point> wave_maker_shape = CreatWaveMakerShape();
		soild_body_part_region_.add_polygon(wave_maker_shape, RegionBooleanOps::add);
		//finish the region modeling
		soild_body_part_region_.done_modeling();

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
	WaveMaking(SolidBody *solid_body, SolidBodyPart *constrained_region)
		: ConstrainSolidBodyRegion(solid_body, constrained_region)
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
	WeaklyCompressibleFluid fluid("Water", water_block, rho0_f, c_f, mu_f, k_f);
	//creat fluid particles
	FluidParticles fluid_particles(water_block);

	//the wall boundary
	WallBoundary *wall_boundary 
		= new WallBoundary(system, "Wall", 0, ParticlesGeneratorOps::lattice);
	//creat solid particles
	SolidParticles solid_particles(wall_boundary);

	//the gate base
	GateBase *gate_base 
		= new GateBase(system, "GateBase", 0, ParticlesGeneratorOps::lattice);
	//elastic soild material properties
	ElasticSolid gate_base_material("ElasticSolid", gate_base, rho0_s, Youngs_modulus, poisson, 0.0);
	//creat particles for the gate base
	ElasticSolidParticles gate_base_particles(gate_base);

	//the elastic gate
	Gate *gate = 
		new Gate(system, "Gate", 0, ParticlesGeneratorOps::lattice);
	//elastic soild material properties
	ElasticSolid gate_material("ElasticSolid", gate, rho0_s, Youngs_modulus, poisson, 0.0);
	//creat particles for the elastic gate
	ElasticSolidParticles gate_particles(gate);

	//set body contact map
	//the contact map gives the data conntections between the bodies
	//basically the the rang of bidies to build neighbor particle lists
	SPHBodyTopology body_topology 
		= { { water_block, { wall_boundary, gate_base, gate } },
		{ wall_boundary, { } },{ gate_base, { gate } },
	    { gate, { gate_base, water_block} } };
	system.SetBodyTopology(&body_topology);

	//setting up the simulation
	system.SetupSPHSimulation();


	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	//methods only used only once
	//-------------------------------------------------------------------

	  /** initial condition for fluid body */
	fluid_dynamics::WeaklyCompressibleFluidInitialCondition set_all_fluid_particles_at_rest(water_block);
	//obtain the initial number density
	fluid_dynamics::InitialNumberDensity
		fluid_initial_number_density(water_block, 
			{ wall_boundary,  gate_base, gate });


	/** initial condition for the solid body */
	solid_dynamics::SolidDynamicsInitialCondition set_all_wall_particles_at_rest(wall_boundary);
	/** initial condition for the elastic solid bodies */
	solid_dynamics::ElasticSolidDynamicsInitialCondition set_all_gate_base_particles_at_rest(gate_base);
	solid_dynamics::ElasticSolidDynamicsInitialCondition set_all_gate_particles_at_rest(gate);
	//initialize normal direction of the wall boundary
	solid_dynamics::NormalDirectionSummation 
		get_wall_normal(wall_boundary, {});
	//initialize normal direction of the gate base
	solid_dynamics::NormalDirectionSummation 
		get_gate_base_normal(gate_base, { gate });
	//initialize normal direction of the elastic gate
	solid_dynamics::NormalDirectionSummation 
		get_gate_normal(gate, { gate_base });

	//corrected strong configuration	
	//gate base
	solid_dynamics::CorrectConfiguration
		gate_base_corrected_configuration_in_strong_form(gate_base, { gate });
	//elastic gate
	solid_dynamics::CorrectConfiguration
		gate_corrected_configuration_in_strong_form(gate, { gate_base });

	//--------------------------------------------------------------------------
	//methods used for time stepping
	//--------------------------------------------------------------------------
	/** Add particle acceleration due to gravity force. */
	//fluid dynamics
	InitializeOtherAccelerations 	initialize_fluid_acceleration(water_block, &gravity);
	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummationFreeSurface
		update_fluid_desnity(water_block, { wall_boundary,  gate_base, gate });
	//time step size without considering sound wave speed
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationVerletFreeSurface
		pressure_relaxation(water_block, 
			{ wall_boundary,  gate_base, gate }, &gravity);

	//FSI
	solid_dynamics::FluidPressureForceOnSolid
		fluid_pressure_force_on_gate(gate, { water_block }, &gravity);

	//solid dynmaics
	//time step size caclutation
	solid_dynamics::GetAcousticTimeStepSize gate_computing_time_step_size(gate);
	//stress relaxation for the gate
	solid_dynamics::StressRelaxation
		gate_stress_relaxation(gate, { gate_base });
	//stress update for contrained wall body
	solid_dynamics::StressInConstrinedElasticBodyFirstHalf
		gate_base_stress_update_first_half(gate_base);
	solid_dynamics::StressInConstrinedElasticBodySecondHalf
		gate_base_stress_update_second_half(gate_base, { gate });

	//constrain region of the inserted body
	WaveMaking
		wave_making(wall_boundary, new WaveMaker(wall_boundary, "WaveMaker"));
	//average velocity
	solid_dynamics::InitializeDisplacement
		gate_initialize_displacement(gate);
	solid_dynamics::UpdateAverageVelocity
		gate_average_velocity(gate);
	solid_dynamics::UpdateElasticNormalDirection 
		gate_update_normal(gate);
	//-------------------------------------------------------------------
	//methods used for updating data structure
	//-------------------------------------------------------------------
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList
		update_water_block_cell_linked_list(water_block);
	//update the configuration of bodies when neccessary
	ParticleDynamicsConfiguration
		update_water_block_configuration(water_block);
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList
		update_wall_boundry_cell_linked_list(wall_boundary);
	//update the configuration of bodies when neccessary
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList
		update_gate_cell_linked_list(gate);
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

	//initial output
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real End_Time = 10.0;
	//time step size for oupt file
	Real D_Time = End_Time/100.0;
	Real Dt = 0.0;//default advection time step sizes
	Real dt = 0.0; //default accoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < D_Time) {

			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			initialize_fluid_acceleration.parallel_exec();
			gate_update_normal.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				wave_making.parallel_exec(dt);
				//fluid dynamics
				pressure_relaxation.parallel_exec(dt);

				//FSI on pressure force
				fluid_pressure_force_on_gate.parallel_exec();

				//solid dynamics
				Real dt_s_sum = 0.0;
				gate_initialize_displacement.parallel_exec();
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
					gate_base_stress_update_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				gate_average_velocity.parallel_exec(dt);

				ite++;
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			
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
