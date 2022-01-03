/**
 * @file 	wetting.h
 * @brief 	Numerical parameters and body definition for 2D two-phase wetting flow.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;						   /**< Tank length. */
Real DH = 1.0;						   /**< Tank height. */
Real particle_spacing_ref = DL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;	   /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;								  /**< Reference density of water. */
Real rho0_a = 1.0e-3;							  /**< Reference density of air. */
Real gravity_g = 0.0;							  /**< Gravity force of fluid. */
Real U_max = 1.0;								  /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;						  /**< Reference sound speed. */
Real mu_f = 5.0e-2;								  /**< Water viscosity. */
Real mu_a = 5.0e-5;								  /**< Air viscosity. */
Real contact_angle = (150.0 / 180.0) * 3.1415926; /**< Contact angle with Wall. */
Real tension_force = 0.008;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
	water_block_shape.push_back(Vecd(0.375 * DL, 0.35 * DH));
	water_block_shape.push_back(Vecd(0.625 * DL, 0.35 * DH));
	water_block_shape.push_back(Vecd(0.625 * DL, 0.0));
	water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
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
/** create inner wall shape */
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
	WaterBlock(SPHSystem &sph_system, const string &body_name)
		: FluidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
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
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
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
		std::vector<Vecd> outer_shape = createOuterWallShape();
		std::vector<Vecd> inner_shape = createInnerWallShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_shape, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
