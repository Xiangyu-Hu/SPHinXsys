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
	Vecd(0.0, 0.0),Vecd(0.0, DH),Vecd(DL, DH),Vecd(DL, 0.0),Vecd(0.0, 0.0)
};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape
{
	Vecd(-BW, -BW),Vecd(-BW, DH + BW),Vecd(DL + BW, DH + BW),Vecd(DL + BW, -BW),Vecd(-BW, -BW),
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape
{
	Vecd(0.0, 0.0),Vecd(0.0, DH),Vecd(DL, DH),Vecd(DL, 0.0),Vecd(0.0, 0.0)
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

class TestVeloGrad: public fluid_dynamics::FluidInitialCondition
{
protected:
	StdLargeVec<Vecd>& pos_, & vel_;
	StdLargeVec<Matd> expect_velo_gradient_;
	StdLargeVec<Matd> &velocity_gradient;
public:
	TestVeloGrad(SPHBody& sph_body): FluidInitialCondition(sph_body),
		pos_(particles_->pos_), vel_(particles_->vel_),
		velocity_gradient(*particles_->getVariableByName<Matd>("Velocity_Gradient"))
	{
		particles_->registerVariable(expect_velo_gradient_,"ExpectVeloGradient");
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
	int getFluidParticleNumber()
	{
		return pos_.size();
	};
	
	Matd getExpectVelocityGradient(size_t index_i)
	{
		return expect_velo_gradient_[index_i];
	};
	Matd getCalculatedVelocityGradient(size_t index_i)
	{
		return velocity_gradient[index_i];
	};
};


TEST(TurbulentModule, TestVelocityGradient)
{
	//MyObject my_obj_(18);
	 //std::cout << "aa" << my_obj_.iii << std::endl;
	BoundingBox system_domain_bounds(Vec2d(-DL - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + BW));
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(sph_system);

	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();

	InnerRelation water_block_inner(water_block);
	ComplexRelation water_block_complex_relation(water_block_inner, { &wall_boundary });

	InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelComplex, SequencedPolicy> k_equation_relaxation(water_block_complex_relation);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_complex_relation);

	BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
	/** Initialize testing. */

	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();

	/** Update cell linked list and configuration. */
	water_block.updateCellLinkedListWithParticleSort(100);
	water_block_complex_relation.updateConfiguration();


	write_body_states.writeToFile();

	update_density_by_summation.exec();
	size_t number_of_iterations=0;
	Real integration_time=0.0;
	Real dt=0.0;
	SimpleDynamics<TestVeloGrad, SequencedPolicy> test_velocity_gradient(water_block);

	test_velocity_gradient.exec(); //** impose initial velocity profile and calculate theoretical velocity gradient value *
	
	k_equation_relaxation.exec(dt);

	write_body_states.writeToFile();


	number_of_iterations++;
	int num_fluid_particle = test_velocity_gradient.getFluidParticleNumber();
	std::cout << "number of fluid particles=" << num_fluid_particle << std::endl;
	for (int index_i = 0; index_i < num_fluid_particle; ++index_i)
	{
		Matd A = test_velocity_gradient.getExpectVelocityGradient(index_i);
		Matd B = test_velocity_gradient.getCalculatedVelocityGradient(index_i);
		for (int j = 0; j < Dimensions; ++j)
		{
			for (int i = 0; i < Dimensions; ++i)
			{
				ASSERT_NEAR(A(i, j), B(i, j), 0.02);
			}
		}
	}
}

//=================================================================================================//
int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}