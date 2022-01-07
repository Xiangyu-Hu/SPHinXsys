/**
 * @file 	two_phase_dambreak.h
 * @brief 	Numerical parameters and body definition for 2D two-phase dambreak flow.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.3;						/**< Tank length. */
Real DH = 2.0;						/**< Tank height. */
Real LL = 2.0;						/**< Liquid colume length. */
Real LH = 1.0;						/**< Liquid colume height. */
Real particle_spacing_ref = 0.05;	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;						 /**< Reference density of water. */
Real rho0_a = 0.001;					 /**< Reference density of air. */
Real gravity_g = 1.0;					 /**< Gravity force of fluid. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.0, 0.0));
	water_block_shape.push_back(Vecd(0.0, LH));
	water_block_shape.push_back(Vecd(LL, LH));
	water_block_shape.push_back(Vecd(LL, 0.0));
	water_block_shape.push_back(Vecd(0.0, 0.0));
	return water_block_shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

	return inner_wall_shape;
}
//----------------------------------------------------------------------
//	Water block body with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &sph_system, const std::string &body_name)
		: FluidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		MultiPolygon original_shape;
		original_shape.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);

		body_shape_.add<MultiPolygonShape>(original_shape);
	}
};
//----------------------------------------------------------------------
//	Air block body with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class AirBlock : public FluidBody
{
public:
	AirBlock(SPHSystem &sph_system, const std::string &body_name)
		: FluidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1.0))
	{
		/** Geomtry definition. */
		MultiPolygon original_shape;
		original_shape.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
		original_shape.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(original_shape);
	}
};
//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		MultiPolygon outer_wall_polygon;
		outer_wall_polygon.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);

		MultiPolygon inner_wall_polygon;
		inner_wall_polygon.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);

		body_shape_.add<MultiPolygonShape>(outer_wall_polygon, "OuterWall");
		body_shape_.substract<MultiPolygonShape>(inner_wall_polygon, "InnerWall");
	}
};
//----------------------------------------------------------------------
//	Observer particle generator.
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(Vecd(DL, 0.2), 0.0));
	}
};
