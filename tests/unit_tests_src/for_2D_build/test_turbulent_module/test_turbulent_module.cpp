#include <gtest/gtest.h>
#include "sphinxsys.h" // SPHinXsys Library.
#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"

using namespace SPH;
Real DL = 0.4;						  /**< Reference length. */
Real DH = 0.2;						  /**< Reference and the height of main channel. */
Real resolution_ref = 0.01;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f;
Real Re = 25000.0;					/**< Reynolds number. */
//Real Re = 100.0;
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */


/** the water block . */
std::vector<Vecd> water_block_shape
{
	Vecd(0.0, 0.0),
	Vecd(0.0, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(0.0, 0.0)
};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape
{
	Vecd(-BW, -BW),
	Vecd(-BW, DH+BW),
	Vecd(DL+BW, DH + BW),
	Vecd(DL + BW, -BW),
	Vecd(-BW, -BW),
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape
{
	Vecd(0.0, 0.0),
	Vecd(0.0, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(0.0, 0.0)
};


class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};

TEST(VelocityGradientCal, VelocityGradient)
{
    Real tolerance = 1e-6;
   







	//EXPECT_NEAR(von_Mises_strain_1, von_Mises_strain_ref_1, tolerance);
    //EXPECT_NEAR(von_Mises_strain_2, von_Mises_strain_ref_2, tolerance);

}
//=================================================================================================//
int main(int argc, char* argv[])
{
	/** Domain bounds of the system. */
	BoundingBox system_domain_bounds(Vec2d(-DL - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + BW));
	SPHSystem system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(system);

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();

	InnerRelation water_block_inner(water_block);
	ComplexRelation water_block_complex_relation(water_block_inner, { &wall_boundary });

	InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelComplex, SequencedPolicy> k_equation_relaxation(water_block_complex_relation);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_complex_relation);
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	
	BodyStatesRecordingToVtp write_body_states(io_environment, system.real_bodies_);

	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_boundary_normal_direction.exec();
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	Real end_time = 100.0;
	Real output_interval = end_time / 20000.0; /**< Time stamps for output of body states. */
	Real dt = 0.0;							 /**< Default acoustic time step sizes. */
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_body_states.writeToFile();
	//----------------------------------------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------------------------------------
	int ITER = 0;
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.exec();
			//Real Dt = get_turbulent_fluid_advection_time_step_size.exec();
			Real Dt = get_fluid_advection_time_step_size.exec();

			update_density_by_summation.exec();

			Real relaxation_time = 0.0;
			dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);

			k_equation_relaxation.exec(dt);

			relaxation_time += dt;
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_complex_relation.updateConfiguration();
		}


		ITER = ITER + 1;
		//std::cout << "ITER=" << ITER << std::endl;
		//if (ITER >=12)
		//{
		//	D_Time = End_Time / 4000.0;
		//	//system("pause");
		//}

	}

	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}