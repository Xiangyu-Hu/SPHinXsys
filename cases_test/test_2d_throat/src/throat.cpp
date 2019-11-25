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
Real gravity_g = 1.0;	/**< Gravity force of fluid. */
Real Re = 1.0;			/**< Reynolds number*/
Real mu_f = sqrt(0.25*rho0_f * powern(0.5*DT, 3)* gravity_g / Re);
Real U_f = 0.25*powern(0.5 * DT, 2)* gravity_g / mu_f;
Real c_f = SMAX(10.0*U_f, 10.0*mu_f/ rho0_f/ particle_spacing_ref);
Real k_f = 0.0;
Real mu_p_f = 0.6*mu_f;
Real lambda_f = 10.0;

//define the fluid body
class FluidBlock : public FluidBody
{
	public:
		FluidBlock(SPHSystem &system, string body_name, 
			int refinement_level, ParticlesGeneratorOps op)
			: FluidBody(system, body_name, refinement_level, op)
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
		int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
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
	Gravity gravity(Vecd(gravity_g, 0.0));

	
	//the water block
	FluidBlock *fluid_block 
		= new FluidBlock(system, "FluidBody", 0, ParticlesGeneratorOps::lattice);
	//fluid material properties
	Oldroyd_B_Fluid fluid("Fluid", fluid_block, rho0_f, c_f, mu_f, k_f, lambda_f, mu_p_f);
	//creat fluid particles
	ViscoelasticFluidParticles fluid_particles(fluid_block);

	//the wall boundary
	WallBoundary *wall_boundary 
		= new WallBoundary(system, "Wall", 0, ParticlesGeneratorOps::lattice);
	//creat solid particles
	SolidParticles solid_particles(wall_boundary);

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
	//initialize normal direction of the wall boundary
	solid_dynamics::NormalDirectionSummation get_wall_normal(wall_boundary, {});

	//-------------------------------------------------------------------
	//methods used for time stepping
	//-------------------------------------------------------------------

	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(fluid_block, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(fluid_block, 0);

	
	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummation
		update_fluid_desnity(fluid_block, { wall_boundary });
	//time step size without considering sound wave speed
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(fluid_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(fluid_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationFirstHalfOldroyd_B
		pressure_relaxation_first_half(fluid_block, { wall_boundary});
	fluid_dynamics::PressureRelaxationSecondHalfOldroyd_B
		pressure_relaxation_second_half(fluid_block, { wall_boundary });

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

	//initial output
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);

	int number_of_iterations = 0;
	int screen_output_interval = 100;
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
				//fluid dynamics
				pressure_relaxation_first_half.parallel_exec(dt);
				pressure_relaxation_second_half.parallel_exec(dt);

				dt = get_fluid_time_step_size.parallel_exec();
				if ((relaxation_time + dt) >= Dt) dt = Dt - relaxation_time;
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

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
