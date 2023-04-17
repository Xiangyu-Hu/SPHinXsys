#include <gtest/gtest.h>
#include "sphinxsys.h" // SPHinXsys Library.
#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"

using namespace SPH;
Real DL = 6;						  /**< Reference length. */
Real DH = 0.2;						  /**< Reference and the height of main channel. */
Real resolution_ref = 0.01;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f;
//Real Re = 25000.0;					/**< Reynolds number. */
Real Re = 100.0;
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */


/** the water block . */
std::vector<Vecd> water_block_shape
{
	Vecd(-DL, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(-DL, 0.0)
};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape
{
	Vecd(-DL - 2.0 * BW, -BW), //1
	Vecd(-DL - 2.0 * BW, DH + BW), //2
	Vecd(DL + 2.0 * BW , DH + BW), //3
	Vecd(DL + 2.0 * BW , -BW), //4
	Vecd(-DL - 2.0 * BW, -BW), //1
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape
{
	Vecd(-DL - 3.0 * BW, 0.0), //1
	Vecd(-DL - 3.0 * BW, DH), //2
	Vecd(DL + 3.0 * BW  , DH), //3
	Vecd(DL + 3.0 * BW , 0.0), //4
	Vecd(-DL - 3.0 * BW, 0.0), //1
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

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();

	InnerRelation water_block_inner(water_block);
	ComplexRelation water_block_complex_relation(water_block_inner, { &wall_boundary });

	InteractionWithUpdate<fluid_dynamics::K_TurtbulentModelComplex, SequencedPolicy> k_equation_relaxation(water_block_complex_relation);

	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}