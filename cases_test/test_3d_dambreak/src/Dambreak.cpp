/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/
/**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
Real particle_spacing_ref = 0.05; //particle spacing
Real BW = particle_spacing_ref * 4; //boundary width

Real DL = 5.366; 						//tank length
Real DH = 2.0; 							//tank height
Real DW = 0.5;							//tank width
Real LL = 2.0; 							//liquid length
Real LH = 1.0; 							//liquid height
Real LW = 0.5; 						//liquid width

//for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0*sqrt(gravity_g * LH);
Real c_f = 10.0*U_f;
Real mu_f = 0.0;
Real k_f = 0.0;

//for initial condition
Real initial_pressure = 0.0;
Vecd intial_velocity(0.0, 0.0, 0.0);
/* resolution which control the quality of polygonalmesh created by geometry system */
int resolution(50);

//define the fluid body
class WaterBlock : public WeaklyCompressibleFluidBody
{
	public:
		WaterBlock(SPHSystem &system, string body_name,
			WeaklyCompressibleFluid* material, 
			WeaklyCompressibleFluidParticles &weakly_compressible_fluid_particles, 
			int refinement_level, ParticlesGeneratorOps op)
			: WeaklyCompressibleFluidBody(system, body_name, material, 
				weakly_compressible_fluid_particles, refinement_level, op)
		{
			Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
			Vecd translation_water = halfsize_water;
			Geometry *geometry_water = new Geometry(halfsize_water, resolution , translation_water);
			body_region_.add_geometry(geometry_water, RegionBooleanOps::add);


			/* the function name, done_modeling, is confused, modify it in the future */
			body_region_.done_modeling();
		}

	void InitialCondition() 
	{
		Particles &base_particles = base_particles_;
		WeaklyCompressibleFluidParticles &fluid_particles
			= weakly_compressible_fluid_particles_;

		for (int i = 0; i < weakly_compressible_fluid_particles_.number_of_particles_; ++i) {
			fluid_particles.fluid_data_[i].p_ = initial_pressure;
			fluid_particles.base_particle_data_[i].vel_n_ = intial_velocity;
			fluid_particles.base_particle_data_[i].dvel_dt_(0);
			fluid_particles.fluid_data_[i].rho_0_
				= material_->ReinitializeRho(initial_pressure);
			fluid_particles.fluid_data_[i].rho_n_
				= material_->ReinitializeRho(initial_pressure);
			fluid_particles.fluid_data_[i].mass_
				= fluid_particles.fluid_data_[i].rho_0_
				*base_particles.base_particle_data_[i].Vol_;
		}
	}
};

//define the static solid wall boudary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, string body_name, 
		SolidBodyParticles &solid_particles, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, solid_particles, refinement_level, op)
	{
		Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
		Vecd translation_wall(0.5 * DL, 0.5 * DH, 0.5 * DW);
		Geometry *geometry_outer = new Geometry(halfsize_outer, resolution, translation_wall);
		body_region_.add_geometry(geometry_outer, RegionBooleanOps::add);

		Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
		Geometry *geometry_inner = new Geometry(halfsize_inner, resolution, translation_wall);
		body_region_.add_geometry(geometry_inner, RegionBooleanOps::sub);

		body_region_.done_modeling();
	}

	void InitialCondition() 
	{
		for (int i = 0; i < solid_particles_.number_of_particles_; ++i) {
			solid_particles_.base_particle_data_[i].vel_n_ = intial_velocity;
			Vecd zero(0);
			solid_particles_.base_particle_data_[i].dvel_dt_ = zero;
			solid_particles_.solid_body_data_[i].vel_ave_ = zero;
			solid_particles_.solid_body_data_[i].dvel_dt_ave_ = zero;
		}
	}
};

//define an observer body
class FluidObserver : public ObserverEulerianBody
{
public:
	FluidObserver(SPHSystem &system, string body_name,
		ObserverParticles &observer_particles, int refinement_level, ParticlesGeneratorOps op)
		: ObserverEulerianBody(system, body_name, observer_particles, refinement_level, op)
	{
		//add observation point
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.01, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.1, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.2, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.24, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.252, 0.5 * DW), 0.0));
		body_input_points_volumes_.push_back(make_pair(Point(DL, 0.266, 0.5 * DW), 0.0));
	}
};

//the main program
int main()
{

	//build up context -- a SPHSystem
	SPHSystem system(Vecd(-BW, -BW, -BW), 
		Vecd(DL + BW, DH + BW, DW + BW), particle_spacing_ref);

	//Configuration of Materials
	WeaklyCompressibleFluid fluid("Water", rho0_f, c_f, mu_f, k_f);
	//creat a fluid particle cotainer
	WeaklyCompressibleFluidParticles fluid_particles("WaterBody");
	//the water block
	WaterBlock *water_block 
		= new WaterBlock(system, "WaterBody", &fluid, fluid_particles, 0, ParticlesGeneratorOps::lattice);
	
	//creat a solid particle cotainer
	SolidBodyParticles solid_particles("Wall");
	//the wall boundary
	WallBoundary *wall_boundary 
		= new WallBoundary(system, "Wall", solid_particles, 0, ParticlesGeneratorOps::lattice);

	//create a observer particle container
	ObserverParticles observer_particles("Fluidobserver");
	FluidObserver *fluid_observer 
		= new FluidObserver(system, "Fluidobserver", observer_particles, 0, ParticlesGeneratorOps::direct);

	//define external force
	Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));


	//set body contact map
	//the contact map gives the data conntections between the bodies
	//basically the the rang of bidies to build neighbor particle lists
	SPHBodyTopology body_topology = { { water_block, { wall_boundary } }, 
		{ wall_boundary, {} },{ fluid_observer,{ water_block} } };
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
	solid_dynamics::NormalDirectionSummation 
		get_wall_normal(wall_boundary, {});
	get_wall_normal.exec();
	//obtain the initial number density
	fluid_dynamics::InitialNumberDensity 
		fluid_initial_number_density(water_block, { wall_boundary });
	fluid_initial_number_density.parallel_exec();


	//-------- common paritcle dynamics ----------------------------------------
	InitializeOtherAccelerations
		initialize_fluid_acceleration(water_block, &gravity);

	//-------- fluid dynamics --------------------------------------------------
	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummationFreeSurface
		update_fluid_desnity(water_block, { wall_boundary });
	//time step size without considering sound wave speed
	fluid_dynamics::GetAdvectionTimeStepSize	get_fluid_adevction_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::GetAcousticTimeStepSize		get_fluid_time_step_size(water_block);

	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationVerletFreeSurface
		pressure_relaxation(water_block, { wall_boundary }, &gravity);

	//--------------------------------------------------------------------------
	//methods used for updating data structure
	//--------------------------------------------------------------------------
	//update the cell linked list of bodies when neccessary
	ParticleDynamicsCellLinkedList
		update_cell_linked_list(water_block);
	//update the configuration of bodies when neccessary
	ParticleDynamicsConfiguration 
		update_particle_configuration(water_block);
	ParticleDynamicsContactConfiguration
		update_observer_contact_configuration(fluid_observer);

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	Output output(system);
	WriteBodyStatesToVtu write_water_block_states(output, system.real_bodies_);
	WriteWaterMechanicalEnergy write_water_mechanical_energy(output, water_block, &gravity);
	WriteWaterFront write_water_front(output, water_block);
	WriteObservedFluidPressure
		write_recorded_water_pressure(output, fluid_observer, { water_block });

	//-------------------------------------------------------------------
	//from here the time stepping begines
	//-------------------------------------------------------------------
	//starting time zero
	GlobalStaticVariables::physical_time_ = 0.0;
	write_water_block_states.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real End_Time = 20.0;
	//time step size for oupt file
	Real D_Time = End_Time/20.0;
	Real Dt = 0.0;//default advection time step sizes
	Real dt = 0.0; //default accoustic time step sizes

	//output for initial particles, global data
	write_water_block_states.WriteToFile(GlobalStaticVariables::physical_time_);
	//output global basic parameters
	output.WriteCaseSetup(End_Time, D_Time, GlobalStaticVariables::physical_time_);

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	
	int count_sorting = 0;
	//computation loop starts 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		//integrate time (loop) until the next output time
		while (integeral_time < D_Time) 
		{

			//acceleration due to viscous force and gravity
			initialize_fluid_acceleration.parallel_exec();
			Dt = get_fluid_adevction_time_step_size.parallel_exec();
			update_fluid_desnity.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{
			
				if (ite % 100 == 0) 
				{
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}

				pressure_relaxation.parallel_exec(dt);

				ite++;
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			
			update_cell_linked_list.parallel_exec();
			update_particle_configuration.parallel_exec();
			update_observer_contact_configuration.parallel_exec();
			write_recorded_water_pressure.WriteToFile(GlobalStaticVariables::physical_time_);

		}

		write_water_mechanical_energy.WriteToFile(GlobalStaticVariables::physical_time_);
		write_water_front.WriteToFile(GlobalStaticVariables::physical_time_);

		tick_count t2 = tick_count::now();
		write_water_block_states.WriteToFile(GlobalStaticVariables::physical_time_ * 0.001);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
