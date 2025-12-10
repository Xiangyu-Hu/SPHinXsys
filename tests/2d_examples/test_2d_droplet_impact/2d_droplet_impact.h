/**
 * @file 	2d_droplet_impact.cpp
 * @brief 	A liquid drop impact the solid boundary.
 * @author Shuaihao Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.

using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real droplet_radius = 0.5;
Real DL = 10*droplet_radius;                         /**< Tank length. */
Real DH = 3*droplet_radius;                         /**< Tank height. */
Real particle_spacing_ref = droplet_radius / 20.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;    /**< Extending width for BCs. */
Vecd droplet_center(DL/2.0, droplet_radius);
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;       /**< Reference density of water. */
Real rho0_a = 1.0e-3;       /**< Reference density of air. */
Real U_ref = 5.0;        /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref; /**< Reference sound speed. */
Real mu_f = 5.0e-3;        /**< Water viscosity. */
Real mu_a = 5.0e-5;        /**< Air viscosity. */
Real surface_tension_coeff = 0.1;
Real Re_water = rho0_f *U_ref*2*droplet_radius/mu_f;
Real We_water = rho0_f * U_ref*U_ref*2*droplet_radius / surface_tension_coeff;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
Vecd air_halfsize = inner_wall_halfsize;
Vecd air_translation = inner_wall_translation;
Vecd droplet_halfsize = Vec2d(droplet_radius, droplet_radius);
Vecd droplet_translation = droplet_center;
//----------------------------------------------------------------------
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBall>(droplet_center, droplet_radius);
    }
};
//----------------------------------------------------------------------
// Air body shape definition.
//----------------------------------------------------------------------/**
class AirBlock : public ComplexShape
{
public:
    explicit AirBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(air_translation), air_halfsize);
        subtract<GeometricShapeBall>(droplet_center, droplet_radius);
    }
};
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
//----------------------------------------------------------------------
//	Initial velocity
//----------------------------------------------------------------------
class InitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
public:
    explicit InitialCondition(RealBody& real_body)
        : fluid_dynamics::FluidInitialCondition(real_body) {};

protected:
    void update(size_t index_i, Real dt)
    {
        vel_[index_i][1] = -U_ref;
    };
};