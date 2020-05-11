/**
 * @file 	fsi.cpp
 * @brief 	This is the benchmark test of fliud-structure interaction.
 * @details We consider a flow-induced vibration of an elastic beam behind a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 * @version 0.1
 */
 /**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"

using namespace SPH;

//for system discretation
Real particle_spacing_ref 	= 0.1; 							//particle spacing
Real BW 				= particle_spacing_ref * 4.0; 		//boundary width
Real DLsponge 			= particle_spacing_ref * 20.0;		//sponge region to impose inflow condition

//for geometry
Real DL = 11.0; 	//channel length
Real DH = 4.0; 		//channel height
Real DW = 4.0; 		//channel width

SimTK::UnitVec3 Cylinder_axis(0.0, 0.0, 1.0);
Real Cylinder_radius = 0.5;	//diameter of cylinder
Real Cylinder_halflength = 1.0;	//length of cylinder
Vecd Cylinder_center(2.0, 0.5 * DH, 0.5 * DW );
int Cylinder_resolution(20);

Vecd Flag_halfsize(4.5 * Cylinder_radius, 0.2 * Cylinder_radius, 1.0);
Vecd Flag_center(2.0 + 4.5 * Cylinder_radius , 0.5 * DH, 0.5 * DW );
int Flag_resolution(20);


//for material properties of the fluid
Real rho0_f 	= 1.0;
Real U_f 		= 1.0;
Real c_f 		= 20.0*U_f;
Real Re 		= 100.0;
Real mu_f 		= rho0_f * U_f * (2.0 * Cylinder_radius) / Re;

//for material properties of the solid
Real rho0_s 	= 10.0; 								//reference density
Real poisson 	= 0.4; 									//Poisson ratio
Real Ae 		= 1.4e3; 								//normalized Youngs Modulus
Real Youngs_modulus 		= Ae * rho0_f * U_f * U_f;

//------------------------------------------------------------------------------
//define geometry and initial conditions of SPH bodies
//------------------------------------------------------------------------------
/**
* @brief define the water block geometry
*/
Geometry *CreateWaterBlock()
{
	Vecd halfsize_water(0.5 * (DLsponge + DL), 0.5 * DH, 0.5 * DW);
	Vecd translation_water(0.5 * (DL - DLsponge), 0.5 * DH, 0.5 * DW);
	Geometry *geometry_water = new Geometry(halfsize_water, 50, translation_water);

	return geometry_water;
}
/**
* @brief define the insert cylinder geometry
*/
Geometry *CreateInsertCylinder()
{
	Geometry *geometry_insert_cylinder = new Geometry(Cylinder_axis, Cylinder_radius,
		Cylinder_halflength, Cylinder_resolution, Cylinder_center);

	return geometry_insert_cylinder;
}
/**
* @brief define the attatched flag geometry
*/
Geometry *CreateAttatchedFlag()
{
	Geometry *geometry_attached_flag 
		= new Geometry(Flag_halfsize, Flag_resolution, Flag_center);

	return geometry_attached_flag;
}
/**
* @brief define the attatched flag geometry
*/
Geometry *CreateInflowBuffer()
{
	Vecd halfsize_inflow(0.5 * DLsponge, 0.5 * DH, 0.5 * DW);
	Vecd translation_inflow(-0.5 * DLsponge, 0.5  * DH, 0.5 * DW);
	Geometry *geometry_infow = new Geometry(halfsize_inflow, 20, translation_inflow);

	return geometry_infow;
}
/**
* @brief define the outer wall geometry
*/
Geometry *CreateOuterWall()
{
	Vecd halfsize_outer(0.5 * (DL + DLsponge) + BW, 0.5 * DH + BW, 0.5 * DW + BW);
	Vecd translation_wall(0.5 * (DL - DLsponge), 0.5 * DH, 0.5 * DW);
	Geometry *geometry_outer = new Geometry(halfsize_outer, 50, translation_wall);

	return geometry_outer;
}
/**
* @brief define the inner wall geometry
*/
Geometry *CreateInnerWall()
{
	Vecd halfsize_inner(0.5 * (DL + DLsponge) + 2.0 * BW, 0.5 * DH, 0.5 * DW);
	Vecd translation_wall(0.5 * (DL - DLsponge), 0.5 * DH, 0.5 * DW);
	Geometry *geometry_inner = new Geometry(halfsize_inner, 50, translation_wall);

	return geometry_inner;
}
//define the fluid body
class WaterBlock : public FluidBody
{
	public:
		WaterBlock(SPHSystem &system, string body_name,
			int refinement_level, ParticlesGeneratorOps op)
			: FluidBody(system, body_name, refinement_level, op)
		{

			body_region_.add_geometry(CreateWaterBlock(), RegionBooleanOps::add);
			body_region_.add_geometry(CreateInsertCylinder(), RegionBooleanOps::sub);
			body_region_.add_geometry(CreateAttatchedFlag(), RegionBooleanOps::sub);
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
		mu_ = mu_f;

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
		body_region_.add_geometry(CreateOuterWall(), RegionBooleanOps::add);
		body_region_.add_geometry(CreateInnerWall(), RegionBooleanOps::sub);
		//finish the region modeling
		body_region_.done_modeling();
	}
};

//insert elastic body with constraint
class InsertedBody : public SolidBody
{
public:
	InsertedBody(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		body_region_.add_geometry(CreateInsertCylinder(), RegionBooleanOps::add);
		body_region_.add_geometry(CreateAttatchedFlag(), RegionBooleanOps::add);
		//finish the region modeling
		body_region_.done_modeling();
	}
};
/**
*@brief Define gate material.
*/
class InsertBodyMaterial : public LinearElasticSolid
{
public:
	InsertBodyMaterial() : LinearElasticSolid()
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
class FlagHolder : public BodyPartByParticle
{
public:
	FlagHolder(SolidBody *solid_body, string constrianed_region_name)
		: BodyPartByParticle(solid_body, constrianed_region_name)
	{
		//geometry
		body_part_region_.add_geometry(CreateInsertCylinder(), RegionBooleanOps::add);
		body_part_region_.add_geometry(CreateAttatchedFlag(), RegionBooleanOps::sub);
		//finish the region modeling
		body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
	}
};
/**
* @brief inflow buffer
*/
class InflowBuffer : public BodyPartByCell
{
public:
	InflowBuffer(FluidBody* fluid_body, string constrianed_region_name)
		: BodyPartByCell(fluid_body, constrianed_region_name)
	{
		/** Geomerty definition. */
		body_part_region_.add_geometry(CreateInflowBuffer(), RegionBooleanOps::add);
		/** Finalize the geometry definition and correspoding opertation. */
		body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartCells();
	}
};
//define an observer body
class Observer : public FictitiousBody
{
public:
	Observer(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{
		//add observation point
		body_input_points_volumes_.push_back(
			make_pair(Point(2.0 + 9.0 * Cylinder_radius, 0.5 * DH, 
				0.5 * DW - Cylinder_halflength), 0.0));
		body_input_points_volumes_.push_back(
			make_pair(Point(2.0 + 9.0 * Cylinder_radius, 
				0.5 * DH, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(
			make_pair(Point(2.0 + 9.0 * Cylinder_radius, 0.5 * DH, 
				0.5 * DW + Cylinder_halflength), 0.0));
	}
};

//inflow boundary condition
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	//parameters for define the inflow profile
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(FluidBody* fluid_body,
		BodyPartByCell*constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region)
	{
		u_ave_ = 0.0;
		u_ref_ = 1.0;
		t_ref = 2.0;
	}

	Vecd GetInflowVelocity(Vecd &position, Vecd &velocity)
	{
		Real u = velocity[0];
		Real v = velocity[1];
		Real w = velocity[2];
		if (position[0] < 0.0) {
			u = 6.0*u_ave_*position[1] * (DH - position[1]) / DH / DH;
			v = 0.0;
			w = 0.0;
		}
		return Vecd(u, v, w);
	}

	void PrepareConstraint() override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref 
			? 0.5*u_ref_*(1.0 - cos(pi*run_time/ t_ref)) : u_ref_;
	}
};

//the main program
int main()
{
	//build up context -- a SPHSystem
	SPHSystem system(Vecd(-DLsponge - BW, -BW, -BW), 
		Vecd(DL + BW, DH + BW, DW + BW), particle_spacing_ref);

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


	//the inseted body immersed in water
	InsertedBody *inserted_body 
		= new InsertedBody(system, "InsertedBody", 1, ParticlesGeneratorOps::lattice);
	//elastic soild material properties
	InsertBodyMaterial *inserted_body_material = new InsertBodyMaterial();
	//creat particles for the elastic body
	ElasticSolidParticles inserted_body_particles(inserted_body, inserted_body_material);

	Observer *flag_observer 
		= new Observer(system, "Observer", 1, ParticlesGeneratorOps::direct);
	//create observer particles 
	BaseParticles observer_particles(flag_observer);

	//set body contact map
	//the contact map gives the data conntections between the bodies
	//basically the the rang of bidies to build neighbor particle lists
	SPHBodyTopology body_topology = { { water_block, { wall_boundary, inserted_body } },
		{ wall_boundary, { } },{ inserted_body, { water_block } }, { flag_observer,{inserted_body} } };
	system.SetBodyTopology(&body_topology);

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	//methods only used only once
	//-------------------------------------------------------------------
	//initialize normal direction of the wall boundary
	solid_dynamics::NormalDirectionSummation get_wall_normal(wall_boundary, {});
	//initialize normal direction of the inserted body
	solid_dynamics::NormalDirectionSummation get_inserted_body_normal(inserted_body, {});
	//corrected strong configuration	
	solid_dynamics::CorrectConfiguration
		inserted_body_corrected_configuration_in_strong_form(inserted_body);

	//-------------------------------------------------------------------
	//methods used for time stepping
	//-------------------------------------------------------------------
	
	//-------- common paritcle dynamics ----------------------------------------
	InitializeATimeStep 	initialize_a_fluid_step(water_block);
	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(water_block, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(water_block, 0);

	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummation
		update_fluid_desnity(water_block, { wall_boundary, inserted_body });
	//time step size without considering sound wave speed
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationFirstHalfRiemann
		pressure_relaxation_first_half(water_block, { wall_boundary, inserted_body });
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(water_block, { wall_boundary, inserted_body });
	//computing viscous acceleration
	fluid_dynamics::ComputingViscousAcceleration
		viscous_acceleration(water_block, { wall_boundary, inserted_body });
	//impose transport velocity
	fluid_dynamics::TransportVelocityCorrection
		transport_velocity_correction(water_block, { wall_boundary, inserted_body });
	//inflow boundary condition
	ParabolicInflow parabolic_inflow(water_block, new InflowBuffer(water_block, "InflowBuffer"));

	//FSI
	solid_dynamics::FluidPressureForceOnSolid
		fluid_pressure_force_on_insrted_body(inserted_body, { water_block });
	solid_dynamics::FluidViscousForceOnSolid
		fluid_viscous_force_on_insrted_body(inserted_body, { water_block });

	//solid dynmaics
	//time step size caclutation
	solid_dynamics::GetAcousticTimeStepSize inserted_body_computing_time_step_size(inserted_body);
	//stress relaxation for the flag
	solid_dynamics::StressRelaxationFirstHalf
		inserted_body_stress_relaxation_first_half(inserted_body);
	solid_dynamics::StressRelaxationSecondHalf
		inserted_body_stress_relaxation_second_half(inserted_body);
	//constrain region of the inserted body
	solid_dynamics::ConstrainSolidBodyRegion
		inserted_body_constrain(inserted_body, new FlagHolder(inserted_body, "FlagHolder"));
	//average velocity
	solid_dynamics::InitializeDisplacement
		inserted_body_initialize_displacement(inserted_body);
	solid_dynamics::UpdateAverageVelocity
		inserted_body_average_velocity(inserted_body);

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
		update_inserted_body_cell_linked_list(inserted_body);
	ParticleDynamicsContactConfiguration
		update_inserted_body_contact_configuration(inserted_body);

	//-------------------------------------------------------------------
	//from here the time stepping begines
	//-------------------------------------------------------------------
	//starting time zero
	GlobalStaticVariables::physical_time_ = 0.0;

	/** Pre-simultion*/

	//initial periodic boundary condition
	//which copies the particle identifies
	//as extra cell linked list form 
	//periodic regions to the corresponding boundaries
	//for buiding up of extra configuration
	system.InitializeSystemCellLinkedLists();
	periodic_condition.parallel_exec();
	//update configuration after periodic boundary condition
	system.InitializeSystemConfigurations();

	get_wall_normal.parallel_exec();
	get_inserted_body_normal.parallel_exec();
	inserted_body_corrected_configuration_in_strong_form.parallel_exec();

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	WriteBodyStatesToVtu write_real_body_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<Vecd, BaseParticles,
		BaseParticleData, &BaseParticles::base_particle_data_, &BaseParticleData::pos_n_>
		write_flag_free_end("Displacement", in_output, flag_observer, inserted_body);

	//initial output
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_flag_free_end.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real End_Time = 200.0;
	//time step size for oupt file
	Real D_Time = 1.0;
	Real Dt = 0.0;//default advection time step sizes
	Real dt = 0.0; //default accoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//computation loop starts 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		//integrate time (loop) until the next output time
		while (integeral_time < D_Time) 
		{
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			//FSI for viscous force
			fluid_viscous_force_on_insrted_body.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				//fluid dynamics
				pressure_relaxation_first_half.parallel_exec(dt);
				//FSI for pressure force
				fluid_pressure_force_on_insrted_body.parallel_exec();
				pressure_relaxation_second_half.parallel_exec(dt);

				//solid dynamics
				Real dt_s_sum = 0.0;
				inserted_body_initialize_displacement.parallel_exec();
				while (dt_s_sum < dt) {

					Real dt_s = inserted_body_computing_time_step_size.parallel_exec();
					if (dt - dt_s_sum < dt_s) dt_s = dt - dt_s_sum;

					if (ite % 100 == 0) {
						cout << "N=" << ite << " Time: "
							<< GlobalStaticVariables::physical_time_ << "	dt_s: "
							<< dt_s << "\n";
					}

					inserted_body_stress_relaxation_first_half.parallel_exec(dt_s);
					inserted_body_constrain.parallel_exec();
					inserted_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				inserted_body_average_velocity.parallel_exec(dt);

				ite++;
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();

			}
			
			//water block confifuration and periodic constion
			periodic_bounding.parallel_exec();
			update_water_block_cell_linked_list.parallel_exec();
			periodic_condition.parallel_exec();
			update_water_block_configuration.parallel_exec();

			//inserted body contact configuration
			update_inserted_body_cell_linked_list.parallel_exec();
			update_inserted_body_contact_configuration.parallel_exec();

			write_flag_free_end.WriteToFile(GlobalStaticVariables::physical_time_);
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
