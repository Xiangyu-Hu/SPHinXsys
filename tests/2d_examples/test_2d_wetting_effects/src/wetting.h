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
Real DL = 2.0;                         /**< Tank length. */
Real DH = 1.0;                         /**< Tank height. */
Real particle_spacing_ref = DL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;    /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                /**< Reference density of water. */
Real rho0_a = 1.0e-3;                             /**< Reference density of air. */
Real U_max = 1.0;                                 /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                          /**< Reference sound speed. */
Real mu_f = 5.0e-2;                               /**< Water viscosity. */
Real mu_a = 5.0e-4;                               /**< Air viscosity. */
Real contact_angle = (150.0 / 180.0) * 3.1415926; /**< Contact angle with Wall. */
Real tension_force = 0.008;
//----------------------------------------------------------------------
//	Geometric elements used in the shape modeling.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
    water_block_shape.push_back(Vecd(0.375 * DL, 0.35 * DH));
    water_block_shape.push_back(Vecd(0.625 * DL, 0.35 * DH));
    water_block_shape.push_back(Vecd(0.625 * DL, 0.0));
    water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
    return water_block_shape;
}

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
//	case-dependent geometric shape of water block.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	case-dependent geometric shape of air block.
//----------------------------------------------------------------------
class AirBlock : public MultiPolygonShape
{
  public:
    explicit AirBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Wall boundary shape definition.
//----------------------------------------------------------------------
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
