/**
 * @file 	droplet.h
 * @brief 	Numerical parameters and body definition for 2D two-phase flow.
 * @author 	Chi Zhang and Xiangyu Hu
 */
/**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 2.0;						   /**< Tank length. */
Real DH = 2.0;						   /**< Tank height. */
Real LL = 1.0;						   /**< Liquid column length. */
Real LH = 1.0;						   /**< Liquid column height. */
Real particle_spacing_ref = DL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 2;	   /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		  /**< Reference density of water. */
Real rho0_a = 0.001;	  /**< Reference density of air. */
Real U_max = 1.0;		  /**< Characteristic velocity. */
Real c_f = 10.0 * U_max; /**< Reference sound speed. */
Real mu_f = 0.2;		  /**< Water viscosity. */
Real mu_a = 0.0002;		  /**< Air viscosity. */
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.5, 0.5));
	water_block_shape.push_back(Vecd(0.5, LH + 0.5));
	water_block_shape.push_back(Vecd(LL + 0.5, LH + 0.5));
	water_block_shape.push_back(Vecd(LL + 0.5, 0.5));
	water_block_shape.push_back(Vecd(0.5, 0.5));
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
/**
*@brief 	Water body shape definition.
*/
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
	}
};
/**
*@brief 	Air body shape definition.
*/
class AirBlock : public MultiPolygonShape
{
public:
	explicit AirBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
	}
};
/**
 * @brief 	Wall boundary shape definition.
 */
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> outer_shape = createOuterWallShape();
		std::vector<Vecd> inner_shape = createInnerWallShape();
		multi_polygon_.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_shape, ShapeBooleanOps::sub);
	}
};
