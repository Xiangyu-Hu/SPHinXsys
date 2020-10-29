/**
 * @file 	fsi.cpp
 * @brief 	This is the benchmark test of fluid-structure interaction.
 * @details We consider a flow-induced vibration of an elastic beam behind a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 * @version 0.1
 */
 /**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"

using namespace SPH;

//global paramteres
Real particle_spacing_ref 	= 0.1; 							//particle spacing
Real BW 				= particle_spacing_ref * 4.0; 		//boundary width
Real DL_sponge 			= particle_spacing_ref * 20.0;		//sponge region to impose inflow condition

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
TriangleMeshShape* CreateWaterBlock()
{
	Vecd halfsize_water(0.5 * (DL_sponge + DL), 0.5 * DH, 0.5 * DW);
	Vecd translation_water(0.5 * (DL - DL_sponge), 0.5 * DH, 0.5 * DW);
	TriangleMeshShape*geometry_water = new TriangleMeshShape(halfsize_water, 50, translation_water);

	return geometry_water;
}
/**
* @brief define the insert cylinder geometry
*/
TriangleMeshShape* CreateInsertCylinder()
{
	TriangleMeshShape*geometry_insert_cylinder 
		= new TriangleMeshShape(Cylinder_axis, Cylinder_radius,
		Cylinder_halflength, Cylinder_resolution, Cylinder_center);

	return geometry_insert_cylinder;
}
/**
* @brief define the attached flag geometry
*/
TriangleMeshShape* CreateAttachedFlag()
{
	TriangleMeshShape* geometry_attached_flag
		= new TriangleMeshShape(Flag_halfsize, Flag_resolution, Flag_center);

	return geometry_attached_flag;
}
/**
* @brief define the attached flag geometry
*/
TriangleMeshShape* CreateInflowBuffer()
{
	Vecd halfsize_inflow(0.5 * DL_sponge, 0.5 * DH, 0.5 * DW);
	Vecd translation_inflow(-0.5 * DL_sponge, 0.5  * DH, 0.5 * DW);
	TriangleMeshShape* geometry_inflow = new TriangleMeshShape(halfsize_inflow, 20, translation_inflow);
	 
	return geometry_inflow;
}
/**
* @brief define the outer wall geometry
*/
TriangleMeshShape* CreateOuterWall()
{
	Vecd halfsize_outer(0.5 * (DL + DL_sponge) + BW, 0.5 * DH + BW, 0.5 * DW + BW);
	Vecd translation_wall(0.5 * (DL - DL_sponge), 0.5 * DH, 0.5 * DW);
	TriangleMeshShape* geometry_outer = new TriangleMeshShape(halfsize_outer, 50, translation_wall);

	return geometry_outer;
}
/**
* @brief define the inner wall geometry
*/
TriangleMeshShape* CreateInnerWall()
{
	Vecd halfsize_inner(0.5 * (DL + DL_sponge) + 2.0 * BW, 0.5 * DH, 0.5 * DW);
	Vecd translation_wall(0.5 * (DL - DL_sponge), 0.5 * DH, 0.5 * DW);
	TriangleMeshShape *geometry_inner = new TriangleMeshShape(halfsize_inner, 50, translation_wall);

	return geometry_inner;
}
//define the fluid body
class WaterBlock : public FluidBody
{
	public:
		WaterBlock(SPHSystem &system, string body_name,	int refinement_level)
			: FluidBody(system, body_name, refinement_level)
		{

			body_shape_ = new ComplexShape(body_name);
			body_shape_->addTriangleMeshShape(CreateWaterBlock(), ShapeBooleanOps::add);
			body_shape_->addTriangleMeshShape(CreateInsertCylinder(), ShapeBooleanOps::sub);
			body_shape_->addTriangleMeshShape(CreateAttachedFlag(), ShapeBooleanOps::sub);
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
	WallBoundary(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addTriangleMeshShape(CreateOuterWall(), ShapeBooleanOps::add);
		body_shape_->addTriangleMeshShape(CreateOuterWall(), ShapeBooleanOps::sub);
	}
};

//insert elastic body with constraint
class InsertedBody : public SolidBody
{
public:
	InsertedBody(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addTriangleMeshShape(CreateInsertCylinder(), ShapeBooleanOps::add);
		body_shape_->addTriangleMeshShape(CreateAttachedFlag(), ShapeBooleanOps::add);
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
	FlagHolder(SolidBody *solid_body, string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		//geometry
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateInsertCylinder(), ShapeBooleanOps::add);
		body_part_shape_->addTriangleMeshShape(CreateAttachedFlag(), ShapeBooleanOps::sub);

		//tag the constrained particle
		tagBodyPart();
	}
};
/**
* @brief inflow buffer
*/
class InflowBuffer : public BodyPartByCell
{
public:
	InflowBuffer(FluidBody* fluid_body, string constrained_region_name)
		: BodyPartByCell(fluid_body, constrained_region_name)
	{
		/** Geometry definition. */
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addTriangleMeshShape(CreateInflowBuffer(), ShapeBooleanOps::add);

		//tag the constrained particle
		tagBodyPart();
	}
};
//define an observer body
class Observer : public FictitiousBody
{
public:
	Observer(SPHSystem &system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
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

	Vecd getTargetVelocity(Vecd &position, Vecd &velocity)
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

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref 
			? 0.5*u_ref_*(1.0 - cos(Pi*run_time/ t_ref)) : u_ref_;
	}
};

//the main program
int main()
{
	//build up context -- a SPHSystem
	SPHSystem system(Vecd(-DL_sponge - BW, -BW, -BW), 
		Vecd(DL + BW, DH + BW, DW + BW), particle_spacing_ref);

	//the water block
	WaterBlock *water_block = new WaterBlock(system, "WaterBody", 0);
	//fluid material properties
	WaterMaterial *water_material = new WaterMaterial();
	//creat fluid particles
	FluidParticles fluid_particles(water_block, water_material);

	//the wall boundary
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall", 0);
	//creat solid particles
	SolidParticles wall_particles(wall_boundary);


	//the inserted body immersed in water
	InsertedBody *inserted_body = new InsertedBody(system, "InsertedBody", 1);
	//elastic solid material properties
	InsertBodyMaterial *inserted_body_material = new InsertBodyMaterial();
	//creat particles for the elastic body
	ElasticSolidParticles inserted_body_particles(inserted_body, inserted_body_material);

	Observer *flag_observer = new Observer(system, "Observer", 1);
	//create observer particles 
	BaseParticles observer_particles(flag_observer);

	/** topology */
	SPHBodyInnerRelation* water_block_inner = new SPHBodyInnerRelation(water_block);
	SPHBodyInnerRelation* inserted_body_inner = new SPHBodyInnerRelation(inserted_body);
	SPHBodyComplexRelation* water_block_complex = new SPHBodyComplexRelation(water_block_inner, { wall_boundary, inserted_body });
	SPHBodyContactRelation* inserted_body_contact = new SPHBodyContactRelation(inserted_body, { water_block });
	SPHBodyContactRelation* flag_observer_contact = new SPHBodyContactRelation(flag_observer, { inserted_body });

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------
	//corrected strong configuration	
	solid_dynamics::CorrectConfiguration
		inserted_body_corrected_configuration(inserted_body_inner);
	//-------------------------------------------------------------------
	//methods used for time stepping
	//-------------------------------------------------------------------
	//-------- common paritcle dynamics ----------------------------------------
	InitializeATimeStep 	initialize_a_fluid_step(water_block);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList 	periodic_condition(water_block, 0);

	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummation
		update_fluid_density(water_block_complex);
	//time step size without considering sound wave speed
	fluid_dynamics::AdvectionTimeStepSize	get_fluid_advection_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(water_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationFirstHalf
		pressure_relaxation_first_half(water_block_complex);
	fluid_dynamics::PressureRelaxationSecondHalfRiemann
		pressure_relaxation_second_half(water_block_complex);
	//computing viscous acceleration
	fluid_dynamics::ViscousAcceleration viscous_acceleration(water_block_complex);
	//impose transport velocity
	fluid_dynamics::TransportVelocityFormulation transport_velocity_formulation(water_block_complex);
	//inflow boundary condition
	ParabolicInflow parabolic_inflow(water_block, new InflowBuffer(water_block, "InflowBuffer"));

	//FSI
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_inserted_body(inserted_body_contact);
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_inserted_body(inserted_body_contact);

	//solid dynamics
	//time step size calculation
	solid_dynamics::AcousticTimeStepSize inserted_body_computing_time_step_size(inserted_body);
	//stress relaxation for the flag
	solid_dynamics::StressRelaxationFirstHalf
		inserted_body_stress_relaxation_first_half(inserted_body_inner);
	solid_dynamics::StressRelaxationSecondHalf
		inserted_body_stress_relaxation_second_half(inserted_body_inner);
	//constrain region of the inserted body
	solid_dynamics::ConstrainSolidBodyRegion
		inserted_body_constrain(inserted_body, new FlagHolder(inserted_body, "FlagHolder"));
	//average velocity
	/** Compute the average velocity on fish body. */
	solid_dynamics::AverageVelocityAndAcceleration	average_velocity_and_acceleration(inserted_body);

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
	//for building up of extra configuration
	system.initializeSystemCellLinkedLists();
	periodic_condition.update_cell_linked_list_.parallel_exec();
	//update configuration after periodic boundary condition
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	inserted_body_particles.initializeNormalDirectionFromGeometry();
	inserted_body_corrected_configuration.parallel_exec();

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	WriteBodyStatesToVtu write_real_body_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::pos_n_>
		write_flag_free_end("Displacement", in_output, flag_observer_contact);

	//initial output
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_flag_free_end.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real End_Time = 200.0;
	//time step size for ouput file
	Real D_Time = 1.0;
	Real Dt = 0.0;//default advection time step sizes
	Real dt = 0.0; //default acoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//computation loop starts 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		//integrate time (loop) until the next output time
		while (integration_time < D_Time) 
		{
			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);

			//FSI for viscous force
			fluid_viscous_force_on_inserted_body.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{

				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
				//fluid dynamics
				pressure_relaxation_first_half.parallel_exec(dt);
				//FSI for pressure force
				fluid_pressure_force_on_inserted_body.parallel_exec();
				pressure_relaxation_second_half.parallel_exec(dt);

				//solid dynamics
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt) 
				{

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
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				ite++;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();

			}
			
			//water block configuration and periodic condition
			periodic_condition.bounding_.parallel_exec();
			water_block->updateCellLinkedList();
			inserted_body->updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
			water_block_complex->updateConfiguration();
			inserted_body_contact->updateConfiguration();
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
