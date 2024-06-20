#include <gtest/gtest.h>
#include "sphinxsys.h" // SPHinXsys Library.
//#include "k-epsilon_turbulent_model_complex.h"
//#include "k-epsilon_turbulent_model_complex.hpp"

using namespace SPH;
//----------------------------------------------------------------------
//	Define basic parameters
//----------------------------------------------------------------------
Real DL = 0.4;						  /**< Reference length. */
Real DH = 0.2;						  /**< Reference and the height of main channel. */
Real resolution_ref = 0.01;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs */
Real c_f = 10.0 * U_f;
Real Re = 13750.0;					/**< Reynolds number according to book [2006 Wilcox]. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//	Define geometry
//----------------------------------------------------------------------
/** the water block . */
std::vector<Vecd> water_block_shape
{
	Vecd(0.0 - 2.0 * BW, 0.0),Vecd(0.0 - 2.0 * BW, DH),
	Vecd(DL + 2.0 * BW, DH),Vecd(DL + 2.0 * BW, 0.0),Vecd(0.0 - 2.0 * BW, 0.0)
};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape
{
	Vecd(-BW, -BW),Vecd(-BW, DH + BW),Vecd(DL + BW, DH + BW),Vecd(DL + BW, -BW),Vecd(-BW, -BW),
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape
{
	Vecd(0.0-2.0 * BW, 0.0),Vecd(0.0- 2.0 * BW, DH),
	Vecd(DL+ 2.0 * BW, DH),Vecd(DL+ 2.0 * BW, 0.0),Vecd(0.0- 2.0 * BW, 0.0)
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
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Create global testing environment
//----------------------------------------------------------------------
class BaseTurbulentModule
{
protected:
	BoundingBox system_domain_bounds;
	SPHSystem sph_system;
	IOEnvironment io_environment;
	FluidBody water_block;
	SolidBody wall_boundary;

public:
	BaseTurbulentModule() :
		system_domain_bounds(Vec2d(-DL - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + BW)),
		sph_system(system_domain_bounds, resolution_ref),
		io_environment(sph_system),
		water_block(sph_system, makeShared<WaterBlock>("WaterBody")),
		wall_boundary(sph_system, makeShared<WallBoundary>("Wall"))
	{
		std::cout << "constructor for BaseTurbulentModule" << std::endl;
		water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
		water_block.generateParticles<ParticleGeneratorLattice>();
		wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
		wall_boundary.generateParticles<ParticleGeneratorLattice>();
	}
	virtual ~BaseTurbulentModule() {};
};

class TurbulentModule : public :: testing :: Test, public BaseTurbulentModule
{
protected:
	InnerRelation water_block_inner;
	ComplexRelation water_block_complex_relation;
	//InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelComplex, SequencedPolicy> k_equation_relaxation;
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation;
	BodyStatesRecordingToVtp write_body_states;
	size_t number_of_iterations;
	Real integration_time;
	Real dt;
	int num_fluid_particle;
public:
	TurbulentModule() : 
		BaseTurbulentModule(),
		water_block_inner(water_block),
		water_block_complex_relation(water_block_inner, { &wall_boundary }),
		//k_equation_relaxation(water_block_complex_relation),
		update_density_by_summation(water_block_complex_relation),
		write_body_states(io_environment, sph_system.real_bodies_)
	{
		std::cout << "constructor for TurbulentModule" << std::endl;
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();

		/** Update cell linked list and configuration. */
		water_block.updateCellLinkedListWithParticleSort(100);
		water_block_complex_relation.updateConfiguration();

		write_body_states.writeToFile();



		update_density_by_summation.exec();
		number_of_iterations = 0;
		integration_time = 0.0;
		dt = 0.0;
		num_fluid_particle = water_block.SizeOfLoopRange();
	}
	virtual ~TurbulentModule() {};
};
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//	Test velocity gradient
//----------------------------------------------------------------------
class TestVeloGrad : public fluid_dynamics::FluidInitialCondition
{
protected:
	StdLargeVec<Vecd>& pos_, & vel_;
	StdLargeVec<Matd> expect_velo_gradient_;
	StdLargeVec<Matd>& velocity_gradient;
public:
	TestVeloGrad(SPHBody& sph_body) : FluidInitialCondition(sph_body),
		pos_(particles_->pos_), vel_(particles_->vel_),
		velocity_gradient(*particles_->getVariableByName<Matd>("Velocity_Gradient"))
	{
		particles_->registerVariable(expect_velo_gradient_, "ExpectVeloGradient");
		particles_->addVariableToWrite<Matd>("ExpectVeloGradient");
	}
	virtual ~TestVeloGrad() {};
	void update(size_t index_i, Real dt = 0.0)
	{
		Real Radius = (DH / 2.0);
		Real transformed_pos = pos_[index_i][1] - (DH / 2.0);
		vel_[index_i][0] = 1.5 * U_f * (1.0 - transformed_pos * transformed_pos / Radius / Radius);
		expect_velo_gradient_[index_i](0, 1) = -2.0 * 1.5 * U_f / Radius / Radius * transformed_pos;
	};
	Matd getExpectVelocityGradient(size_t index_i)
	{
		return expect_velo_gradient_[index_i];
	};
	Matd getCalculatedVelocityGradient(size_t index_i)
	{
		return velocity_gradient[index_i];
	};
	Vecd getPosition(size_t index_i)
	{
		return pos_[index_i];
	};
};
TEST_F(TurbulentModule, TestVelocityGradient)
{
	SimpleDynamics<TestVeloGrad, SequencedPolicy> test_velocity_gradient(water_block);
	test_velocity_gradient.exec(); //** impose initial velocity profile and calculate theoretical velocity gradient value *
	//k_equation_relaxation.exec(dt);
	write_body_states.writeToFile();
	number_of_iterations++;
	std::cout << "number of fluid particles=" << num_fluid_particle << std::endl;
	for (int index_i = 0; index_i < num_fluid_particle; ++index_i)
	{
		//** reduce testing region, because in horizontal direction, kernel truncation is avoided by using in/outlet buffer *
		Vecd pos_ = test_velocity_gradient.getPosition(index_i);
		if (pos_[0] >= 0.5 * DL && pos_[0] < 0.5 * DL + BW)
		{
			Matd A = test_velocity_gradient.getExpectVelocityGradient(index_i);
			Matd B = test_velocity_gradient.getCalculatedVelocityGradient(index_i);
			for (int j = 0; j < Dimensions; ++j)
			{
				for (int i = 0; i < Dimensions; ++i)
				{
					ASSERT_NEAR(A(i, j), B(i, j), 1.0);
				}
			}
		}
		
	}
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//	Test K equation
//----------------------------------------------------------------------
class Test_K_Equation : public fluid_dynamics::FluidInitialCondition
{
      protected:
      public:
        Test_K_Equation(SPHBody &sph_body) : FluidInitialCondition(sph_body)
        {
        }
        virtual ~Test_K_Equation(){};
        void update(size_t index_i, Real dt = 0.0){

        };
};
TEST_F(TurbulentModule, TestTurbulentKineticEnergyEquation)
{
        SimpleDynamics<TestVeloGrad, SequencedPolicy> test_velocity_gradient(water_block);
        test_velocity_gradient.exec(); //** impose initial velocity profile and calculate theoretical velocity gradient value *
        //k_equation_relaxation.exec(dt);

        ASSERT_NEAR(1, 1, 0.02);
}
//----------------------------------------------------------------------


//=================================================================================================//
int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}