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
Real mu_p_f = 0.6 * mu_f;
Real lambda_f = 10.0;

//define the fluid body
class FluidBlock : public FluidBody
{
	public:
		FluidBlock(SPHSystem &system, string body_name, int refinement_level)
			: FluidBody(system, body_name, refinement_level)
		{
			std::vector<Point> pnts;
			pnts.push_back(Point(-0.5*DL, -0.5*DH));
			pnts.push_back(Point(-0.5*DL, 0.5*DH));
			pnts.push_back(Point(-DL / 6.0, 0.5*DH));
			pnts.push_back(Point(-DL / 6.0, - 0.5*DH));
			pnts.push_back(Point(-0.5*DL, -0.5*DH));

			std::vector<Point> pnts1;
			pnts1.push_back(Point(-DL/6.0 - BW, -0.5*DT));
			pnts1.push_back(Point(-DL / 6.0 - BW, 0.5*DT));
			pnts1.push_back(Point(DL / 6.0 + BW, 0.5*DT));
			pnts1.push_back(Point(DL / 6.0 + BW, - 0.5*DT));
			pnts1.push_back(Point(-DL / 6.0 - BW, -0.5*DT));

			std::vector<Point> pnts2;
			pnts2.push_back(Point(DL/6.0, -0.5*DH));
			pnts2.push_back(Point(DL / 6.0, 0.5*DH));
			pnts2.push_back(Point(0.5*DL, 0.5*DH));
			pnts2.push_back(Point(0.5*DL, -0.5*DH));
			pnts2.push_back(Point(DL / 6.0, -0.5*DH));

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
		rho_0_ = rho0_f;
		c_0_ = c_f;
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
	WallBoundary(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		std::vector<Point> pnts3;
		pnts3.push_back(Point(-0.5*DL - BW, -0.5*DH - BW));
		pnts3.push_back(Point(-0.5*DL - BW, 0.5*DH + BW));
		pnts3.push_back(Point(0.5*DL + BW, 0.5*DH + BW));
		pnts3.push_back(Point(0.5*DL + BW, -0.5*DH - BW));
		pnts3.push_back(Point(-0.5*DL - BW, -0.5*DH - BW));

		std::vector<Point> pnts;
		pnts.push_back(Point(-0.5*DL - 2.0*BW, -0.5*DH));
		pnts.push_back(Point(-0.5*DL - 2.0*BW, 0.5*DH));
		pnts.push_back(Point(-DL / 6.0, 0.5*DH));
		pnts.push_back(Point(-DL / 6.0, -0.5*DH));
		pnts.push_back(Point(-0.5*DL - 2.0*BW, -0.5*DH));

		std::vector<Point> pnts1;
		pnts1.push_back(Point(-DL / 6.0 - BW, -0.5*DT));
		pnts1.push_back(Point(-DL / 6.0 - BW, 0.5*DT));
		pnts1.push_back(Point(DL / 6.0 + BW, 0.5*DT));
		pnts1.push_back(Point(DL / 6.0 + BW, -0.5*DT));
		pnts1.push_back(Point(-DL / 6.0 - BW, -0.5*DT));

		std::vector<Point> pnts2;
		pnts2.push_back(Point(DL / 6.0, -0.5*DH));
		pnts2.push_back(Point(DL / 6.0, 0.5*DH));
		pnts2.push_back(Point(0.5*DL + 2.0*BW, 0.5*DH));
		pnts2.push_back(Point(0.5*DL + 2.0*BW, -0.5*DH));
		pnts2.push_back(Point(DL / 6.0, -0.5*DH));

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
	SPHSystem system(Vec2d(-0.5*DL - BW, -0.5*DH - BW), 
		Vec2d(0.5*DL + BW, 0.5*DH + BW), particle_spacing_ref);

	//define external force
	Gravity gravity(Vecd(gravity_g, 0.0));

	
	//the water block
	FluidBlock *fluid_block = new FluidBlock(system, "FluidBody", 0);
	//fluid material properties
	NonNewtonianMaterial *non_newtonian_material = new NonNewtonianMaterial();
	//creat fluid particles
	ViscoelasticFluidParticles fluid_particles(fluid_block, non_newtonian_material);

	//the wall boundary
	WallBoundary *wall_boundary = new WallBoundary(system, "Wall", 0);
	//creat solid particles
	SolidParticles solid_particles(wall_boundary);
	/** topology */
	SPHBodyInnerRelation* fluid_block_inner = new SPHBodyInnerRelation(fluid_block);
	SPHBodyComplexRelation* fluid_block_complex = new SPHBodyComplexRelation(fluid_block_inner, { wall_boundary });
	SPHBodyComplexRelation* wall_complex = new SPHBodyComplexRelation(wall_boundary, {});

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	//methods only used only once
	//-------------------------------------------------------------------
	//initialize normal direction of the wall boundary
	solid_dynamics::NormalDirectionSummation get_wall_normal(wall_complex);

	//-------------------------------------------------------------------
	//methods used for time stepping
	//-------------------------------------------------------------------

	/** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(fluid_block, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(fluid_block, 0);

	
	//evaluation of density by summation approach
	fluid_dynamics::DensityBySummation
		update_fluid_density(fluid_block_complex);
	//time step size without considering sound wave speed
	fluid_dynamics::AdvectionTimeStepSize	get_fluid_advection_time_step_size(fluid_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(fluid_block);
	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationFirstHalfOldroyd_B
		pressure_relaxation_first_half(fluid_block_complex);
	fluid_dynamics::PressureRelaxationSecondHalfOldroyd_B
		pressure_relaxation_second_half(fluid_block_complex);

	//-------- common particle dynamics ----------------------------------------
	InitializeATimeStep 	initialize_a_fluid_step(fluid_block, &gravity);
	//computing viscous acceleration
	fluid_dynamics::ViscousAcceleration viscous_acceleration(fluid_block_complex);
	//impose transport velocity
	fluid_dynamics::TransportVelocityFormulation transport_velocity_formulation(fluid_block_complex);
	//computing vorticity in the flow
	fluid_dynamics::VorticityInFluidField
		compute_vorticity(fluid_block_inner);
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
	//for building up of extra configuration
	system.initializeSystemCellLinkedLists();
	periodic_condition.parallel_exec();
	system.initializeSystemConfigurations();

	//prepare quantities will be used once only
	get_wall_normal.parallel_exec();

	//initial output
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);

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
		while (integration_time < D_Time) {

			initialize_a_fluid_step.parallel_exec();
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_fluid_density.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_formulation.correction_.parallel_exec(Dt);

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				//fluid dynamics
				pressure_relaxation_first_half.parallel_exec(dt);
				pressure_relaxation_second_half.parallel_exec(dt);

				dt = get_fluid_time_step_size.parallel_exec();
				if ((relaxation_time + dt) >= Dt) dt = Dt - relaxation_time;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			//water block configuration and periodic condition
			periodic_bounding.parallel_exec();
			system.initializeSystemCellLinkedLists();
			periodic_condition.parallel_exec();
			system.initializeSystemConfigurations();
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
