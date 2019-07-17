/**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
Real DH = 4.0; //channel height
Real DT = 1.0; //throat height
Real DL = 24.0; //channel length
Real particle_spacing_ref = 0.1; //particle spacing
Real BW = particle_spacing_ref * 4.0; //boundary width

//for material properties of the fluid
Real rho0_f = 1.0;
Real U_f = 1.0;
Real c_f = 10.0*U_f;
Real Re = 10.0;
Real mu_f = rho0_f * U_f * (2.0 * DT) / Re;
Real k_f = 0.0;
Real mu_p_f = 0.6*mu_f;
Real lambda_f = 10.0;

//for fluid initial condition
Real initial_pressure = 0.0;
Vec2d intial_velocity(0.0, 0.0);

//define the fluid body
class FluidBlock : public FluidBody
{
	public:
		FluidBlock(SPHSystem &system, string body_name, Oldroyd_B_Fluid &oldroyd_b_material,
			ViscoelasticFluidParticles &viscoelastic_fluid_particles,
									int refinement_level, ParticlesGeneratorOps op)
			: FluidBody(system, body_name, oldroyd_b_material,
				viscoelastic_fluid_particles, refinement_level, op)
		{
			std::vector<Point> pnts;
			pnts.push_back(Point(-0.5*DL, -0.5*DH));
			pnts.push_back(Point(-0.5*DL, 0.5*DH));
			pnts.push_back(Point(-DL / 6.0, 0.5*DH));
			pnts.push_back(Point(-DL / 6.0, - 0.5*DH));
			pnts.push_back(Point(-0.5*DL, -0.5*DH));
			body_region_.add_polygon(pnts, RegionBooleanOps::add);

			std::vector<Point> pnts1;
			pnts1.push_back(Point(-DL/6.0 - BW, -0.5*DT));
			pnts1.push_back(Point(-DL / 6.0 - BW, 0.5*DT));
			pnts1.push_back(Point(DL / 6.0 + BW, 0.5*DT));
			pnts1.push_back(Point(DL / 6.0 + BW, - 0.5*DT));
			pnts1.push_back(Point(-DL / 6.0 - BW, -0.5*DT));
			body_region_.add_polygon(pnts1, RegionBooleanOps::add);

			std::vector<Point> pnts2;
			pnts2.push_back(Point(DL/6.0, -0.5*DH));
			pnts2.push_back(Point(DL / 6.0, 0.5*DH));
			pnts2.push_back(Point(0.5*DL, 0.5*DH));
			pnts2.push_back(Point(0.5*DL, -0.5*DH));
			pnts2.push_back(Point(DL / 6.0, -0.5*DH));
			body_region_.add_polygon(pnts2, RegionBooleanOps::add);

			//finish the region modeling
			body_region_.done_modeling();
		}
};

//define the static solid wall boudary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name, 
		SolidParticles &solid_particles, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, *(new Solid("EmptyWallMaterial")), solid_particles,	refinement_level, op)
	{
		std::vector<Point> pnts3;
		pnts3.push_back(Point(-0.5*DL - BW, -0.5*DH - BW));
		pnts3.push_back(Point(-0.5*DL - BW, 0.5*DH + BW));
		pnts3.push_back(Point(0.5*DL + BW, 0.5*DH + BW));
		pnts3.push_back(Point(0.5*DL + BW, -0.5*DH - BW));
		pnts3.push_back(Point(-0.5*DL - BW, -0.5*DH - BW));
		body_region_.add_polygon(pnts3, RegionBooleanOps::add);

		std::vector<Point> pnts;
		pnts.push_back(Point(-0.5*DL - 2.0*BW, -0.5*DH));
		pnts.push_back(Point(-0.5*DL - 2.0*BW, 0.5*DH));
		pnts.push_back(Point(-DL / 6.0, 0.5*DH));
		pnts.push_back(Point(-DL / 6.0, -0.5*DH));
		pnts.push_back(Point(-0.5*DL - 2.0*BW, -0.5*DH));
		body_region_.add_polygon(pnts, RegionBooleanOps::sub);

		std::vector<Point> pnts1;
		pnts1.push_back(Point(-DL / 6.0 - BW, -0.5*DT));
		pnts1.push_back(Point(-DL / 6.0 - BW, 0.5*DT));
		pnts1.push_back(Point(DL / 6.0 + BW, 0.5*DT));
		pnts1.push_back(Point(DL / 6.0 + BW, -0.5*DT));
		pnts1.push_back(Point(-DL / 6.0 - BW, -0.5*DT));
		body_region_.add_polygon(pnts1, RegionBooleanOps::sub);

		std::vector<Point> pnts2;
		pnts2.push_back(Point(DL / 6.0, -0.5*DH));
		pnts2.push_back(Point(DL / 6.0, 0.5*DH));
		pnts2.push_back(Point(0.5*DL + 2.0*BW, 0.5*DH));
		pnts2.push_back(Point(0.5*DL + 2.0*BW, -0.5*DH));
		pnts2.push_back(Point(DL / 6.0, -0.5*DH));
		body_region_.add_polygon(pnts2, RegionBooleanOps::sub);

		//finish the region modeling
		body_region_.done_modeling();

	}
};

//the main program
int main()
{
	//build up context -- a SPHSystem
	SPHSystem system(Vec2d(-0.5*DL - BW, -0.5*DH - BW), 
		Vec2d(0.5*DL + BW, 0.5*DH + BW), particle_spacing_ref);

	//define external force
	Gravity gravity(Vecd(1.0, 0.0));

	//fluid material properties
	Oldroyd_B_Fluid fluid("Fluid", rho0_f, c_f, mu_f, k_f, lambda_f, mu_p_f);
	
	//creat a fluid particle cotainer
	ViscoelasticFluidParticles fluid_particles("FluidBody");
	//the water block
	FluidBlock *fluid_block 
		= new FluidBlock(system, "FluidBody", fluid, fluid_particles, 0, ParticlesGeneratorOps::lattice);
	
	//creat a solid particle cotainer
	SolidParticles solid_particles("Wall");
	//the wall boundary
	WallBoundary *wall_boundary 
		= new WallBoundary(system, "Wall", solid_particles, 0, ParticlesGeneratorOps::lattice);

	//set body contact map
	//the contact map gives the data conntections between the bodies
	//basically the the rang of bidies to build neighbor particle lists
	SPHBodyTopology body_topology 
		= { { fluid_block, { wall_boundary} },	{ wall_boundary, { } } };
	system.SetBodyTopology(&body_topology);

	//setting up the simulation
	system.SetupSPHSimulation();

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	//methods only used only once
	//-------------------------------------------------------------------
	/** initial condition */
	fluid_dynamics::Oldroyd_B_FluidInitialCondition set_all_fluid_particles_at_rest(fluid_block);
	//obtain the initial number density
	fluid_dynamics::InitialNumberDensity
		fluid_initial_number_density(fluid_block, {wall_boundary});


	/** initial condition for the solid body */
	solid_dynamics::SolidDynamicsInitialCondition set_all_wall_particles_at_rest(wall_boundary);
	//initialize normal direction of the wall boundary
	solid_dynamics::NormalDirectionSummation get_wall_normal(wall_boundary, {});

	//-------------------------------------------------------------------
	//methods used for time stepping
	//-------------------------------------------------------------------

	//periodic bounding
	PeriodicBoundingInXDirection
		periodic_bounding(fluid_block);
	//periodic boundary condition
	PeriodicConditionInXDirection
		periodic_condition(fluid_block);
	
	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummation
		update_fluid_desnity(fluid_block, { wall_boundary });
	//time step size without considering sound wave speed
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(fluid_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(fluid_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::VerletOldroyd_B_Fluid
		pressure_relaxation(fluid_block, { wall_boundary}, &gravity);

	//-------- common paritcle dynamics ----------------------------------------
	InitializeOtherAccelerations
		initialize_other_acceleration(fluid_block, &gravity);
	//computing viscous acceleration
	fluid_dynamics::ComputingViscousAcceleration
		viscous_acceleration(fluid_block, { wall_boundary});
	//impose transport velocity
	fluid_dynamics::TransportVelocityCorrection
		transport_velocity_correction(fluid_block, { wall_boundary });
	//computing vorticity in the flow
	fluid_dynamics::ComputingVorticityInFluidField
		compute_vorticity(fluid_block);

	//-------------------------------------------------------------------
	//methods used for updating data structure
	//-------------------------------------------------------------------
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList
		update_water_block_cell_linked_list(fluid_block);
	//update the configuration of bodies when neccessary
	ParticleDynamicsConfiguration
		update_water_block_configuration(fluid_block);

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
	
	/** apply initial condition */
	set_all_fluid_particles_at_rest.exec();
	set_all_wall_particles_at_rest.exec();

	//initial periodic boundary condition
	//which copies the particle identifies
	//as extra cell linked list form 
	//periodic regions to the corresponding boundaries
	//for buiding up of extra configuration
	periodic_condition.parallel_exec();
	//update configuration after periodic boundary condition
	update_water_block_configuration.parallel_exec();


	//prepare quantities will be used once only
	get_wall_normal.parallel_exec();
	fluid_initial_number_density.parallel_exec();

	//initial output
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real End_Time = 200.0*5.0;
	//time step size for oupt file
	Real D_Time = End_Time/200.0/5.0;
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
		while (integeral_time < D_Time) {

			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();
			initialize_other_acceleration.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 50 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				//fluid dynamics
				pressure_relaxation.parallel_exec(dt);

				ite++;
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

			}
			
			//water block confifuration and periodic constion
			periodic_bounding.parallel_exec();
			update_water_block_cell_linked_list.parallel_exec();
			periodic_condition.parallel_exec();
			update_water_block_configuration.parallel_exec();
		}

		tick_count t2 = tick_count::now();
		compute_vorticity.parallel_exec();
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
