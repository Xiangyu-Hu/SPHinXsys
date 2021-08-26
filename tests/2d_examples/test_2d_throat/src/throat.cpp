/**
 * @file 	throat.cpp
 * @brief 	2D in a channel with a throat.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for non-Newtonian low Reynolds number flows.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
   */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
Real DH = 4.0; //channel height
Real DT = 1.0; //throat height
Real DL = 24.0; //channel length
Real resolution_ref = 0.1; //particle spacing
Real BW = resolution_ref * 4.0; //boundary width
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-0.5 * DL - BW, -0.5 * DH - BW),
	Vec2d(0.5 * DL + BW, 0.5 * DH + BW));

//for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;	/**< Gravity force of fluid. */
Real Re = 0.1;			/**< Reynolds number*/
Real mu_f = sqrt(0.25*rho0_f * powerN(0.5*DT, 3)* gravity_g / Re);
Real U_f = 0.25*powerN(0.5 * DT, 2)* gravity_g / mu_f;
// For low Reynolds number flow the weakly compressible formulation need to 
// consider viscosity for artificial sound speed.
Real c_f = 10.0 * SMAX(U_f, 2.0 * mu_f / rho0_f / DT);
Real mu_p_f = 0.6 * mu_f;
Real lambda_f = 10.0;

//define the fluid body
class FluidBlock : public FluidBody
{
	public:
		FluidBlock(SPHSystem &system, std::string body_name)
			: FluidBody(system, body_name)
		{
			std::vector<Vecd> pnts;
			pnts.push_back(Vecd(-0.5*DL, -0.5*DH));
			pnts.push_back(Vecd(-0.5*DL, 0.5*DH));
			pnts.push_back(Vecd(-DL / 6.0, 0.5*DH));
			pnts.push_back(Vecd(-DL / 6.0, - 0.5*DH));
			pnts.push_back(Vecd(-0.5*DL, -0.5*DH));

			std::vector<Vecd> pnts1;
			pnts1.push_back(Vecd(-DL/6.0 - BW, -0.5*DT));
			pnts1.push_back(Vecd(-DL / 6.0 - BW, 0.5*DT));
			pnts1.push_back(Vecd(DL / 6.0 + BW, 0.5*DT));
			pnts1.push_back(Vecd(DL / 6.0 + BW, - 0.5*DT));
			pnts1.push_back(Vecd(-DL / 6.0 - BW, -0.5*DT));

			std::vector<Vecd> pnts2;
			pnts2.push_back(Vecd(DL/6.0, -0.5*DH));
			pnts2.push_back(Vecd(DL / 6.0, 0.5*DH));
			pnts2.push_back(Vecd(0.5*DL, 0.5*DH));
			pnts2.push_back(Vecd(0.5*DL, -0.5*DH));
			pnts2.push_back(Vecd(DL / 6.0, -0.5*DH));

			body_shape_ = new ComplexShape(body_name);
			body_shape_->addAPolygon(pnts, ShapeBooleanOps::add);
			body_shape_->addAPolygon(pnts1, ShapeBooleanOps::add);
			body_shape_->addAPolygon(pnts2, ShapeBooleanOps::add);
		}
};
/**
 * @brief 	Case dependent material properties definition.
 */
class NonNewtonianMaterial : public Oldroyd_B_Fluid
{
public:
	NonNewtonianMaterial() : Oldroyd_B_Fluid()
	{
		rho0_ = rho0_f;
		c0_ = c_f;
		mu_ = mu_f;
		mu_p_ =  mu_p_f;
		lambda_ = lambda_f;

		assignDerivedMaterialParameters();
	}
};
//define the static solid wall boundary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		std::vector<Vecd> pnts3;
		pnts3.push_back(Vecd(-0.5*DL - BW, -0.5*DH - BW));
		pnts3.push_back(Vecd(-0.5*DL - BW, 0.5*DH + BW));
		pnts3.push_back(Vecd(0.5*DL + BW, 0.5*DH + BW));
		pnts3.push_back(Vecd(0.5*DL + BW, -0.5*DH - BW));
		pnts3.push_back(Vecd(-0.5*DL - BW, -0.5*DH - BW));

		std::vector<Vecd> pnts;
		pnts.push_back(Vecd(-0.5*DL - 2.0*BW, -0.5*DH));
		pnts.push_back(Vecd(-0.5*DL - 2.0*BW, 0.5*DH));
		pnts.push_back(Vecd(-DL / 6.0, 0.5*DH));
		pnts.push_back(Vecd(-DL / 6.0, -0.5*DH));
		pnts.push_back(Vecd(-0.5*DL - 2.0*BW, -0.5*DH));

		std::vector<Vecd> pnts1;
		pnts1.push_back(Vecd(-DL / 6.0 - BW, -0.5*DT));
		pnts1.push_back(Vecd(-DL / 6.0 - BW, 0.5*DT));
		pnts1.push_back(Vecd(DL / 6.0 + BW, 0.5*DT));
		pnts1.push_back(Vecd(DL / 6.0 + BW, -0.5*DT));
		pnts1.push_back(Vecd(-DL / 6.0 - BW, -0.5*DT));

		std::vector<Vecd> pnts2;
		pnts2.push_back(Vecd(DL / 6.0, -0.5*DH));
		pnts2.push_back(Vecd(DL / 6.0, 0.5*DH));
		pnts2.push_back(Vecd(0.5*DL + 2.0*BW, 0.5*DH));
		pnts2.push_back(Vecd(0.5*DL + 2.0*BW, -0.5*DH));
		pnts2.push_back(Vecd(DL / 6.0, -0.5*DH));

		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(pnts3, ShapeBooleanOps::add);
		body_shape_->addAPolygon(pnts, ShapeBooleanOps::sub);
		body_shape_->addAPolygon(pnts1, ShapeBooleanOps::sub);
		body_shape_->addAPolygon(pnts2, ShapeBooleanOps::sub);
	}
};

//the main program
int main()
{
	//build up context -- a SPHSystem
	SPHSystem system(system_domain_bounds, resolution_ref);

	//define external force
	Gravity gravity(Vecd(gravity_g, 0.0));

	
	//the water block
	FluidBlock *fluid_block = new FluidBlock(system, "FluidBody");
	//fluid material properties
	NonNewtonianMaterial *non_newtonian_material = new NonNewtonianMaterial();
	//creat fluid particles
	ViscoelasticFluidParticles fluid_particles(fluid_block, non_newtonian_material);

	//the wall boundary
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall");
	//creat solid particles
	SolidParticles wall_particles(wall_boundary);
	/** topology */
	BodyRelationInner* fluid_block_inner = new BodyRelationInner(fluid_block);
	ComplexBodyRelation* fluid_block_complex = new ComplexBodyRelation(fluid_block_inner, { wall_boundary });
	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------
	//-------------------------------------------------------------------
	//methods used for time stepping
	//-------------------------------------------------------------------
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirectionUsingGhostParticles 	periodic_condition(fluid_block, xAxis);

	
	//evaluation of density by summation approach
	fluid_dynamics::DensitySummationComplex	update_density_by_summation(fluid_block_complex);
	//time step size without considering sound wave speed and viscosity
	fluid_dynamics::AdvectionTimeStepSizeForImplicitViscosity	get_fluid_advection_time_step_size(fluid_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(fluid_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationRiemannWithWallOldroyd_B	pressure_relaxation(fluid_block_complex);
	pressure_relaxation.pre_processes_.push_back(&periodic_condition.ghost_update_);
	fluid_dynamics::DensityRelaxationRiemannWithWallOldroyd_B	density_relaxation(fluid_block_complex);
	density_relaxation.pre_processes_.push_back(&periodic_condition.ghost_update_);

	//-------- common particle dynamics ----------------------------------------
	TimeStepInitialization 	initialize_a_fluid_step(fluid_block, &gravity);
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(fluid_block_complex);
	//computing viscous effect implicitly and with update velocity directly other than viscous acceleration
	DampingPairwiseWithWall<indexVector, Vec2d, DampingPairwiseInner>
		implicit_viscous_damping(fluid_block_complex, "Velocity", mu_f);

	//impose transport velocity
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(fluid_block_complex);
	//computing vorticity in the flow
	fluid_dynamics::VorticityInner compute_vorticity(fluid_block_inner);
	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_real_body_states(in_output, system.real_bodies_);

	//-------------------------------------------------------------------
	//from here the time stepping begines
	//-------------------------------------------------------------------
	//starting time zero
	GlobalStaticVariables::physical_time_ = 0.0;
	
	system.initializeSystemCellLinkedLists();
	//initial periodic boundary condition
	periodic_condition.ghost_creation_.parallel_exec();
	system.initializeSystemConfigurations();

	//prepare quantities will be used once only
	wall_particles.initializeNormalDirectionFromGeometry();

	//initial output
	write_real_body_states.writeToFile(0);

	int number_of_iterations = 0;
	int screen_output_interval = 100;
	Real End_Time = 100.0;
	//time step size for ouput file
	Real D_Time = End_Time/100.0;
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
			update_density_by_summation.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt);
				pressure_relaxation.parallel_exec(dt);
				implicit_viscous_damping.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			//water block configuration and periodic condition
			periodic_condition.bounding_.parallel_exec();
			system.initializeSystemCellLinkedLists();
			periodic_condition.ghost_creation_.parallel_exec();
			system.initializeSystemConfigurations();
		}

		tick_count t2 = tick_count::now();
		compute_vorticity.parallel_exec();
		write_real_body_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
