/**
 * @file 	2d_eulerian_flow_around_cylinder.h
 * @brief 	This is a test to show the standard eulerian flow around cylinder case.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#ifndef EULERIAN_FLOW_AROUND_CYLINDER_H
#define EULERIAN_FLOW_AROUND_CYLINDER_H
#include "common_shared_eulerian_classes.h"                // shared eulerian classes for weakly-compressible and compressible fluid.
#include "common_weakly_compressible_eulerian_classes.hpp" // eulerian classes for weakly compressible fluid only.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 15.0;                        /**< Channel length. */
Real DH = 10.0;                        /**< Channel height. */
Real resolution_ref = 1.0 / 4.0;       /**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Vec2d cylinder_center(4.0, DH / 2.0);  /**< Location of the cylinder center. */
Real cylinder_radius = 1.0;            /**< Radius of the cylinder. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon circle(cylinder_center, cylinder_radius, 100);
        subtract<MultiPolygonShape>(circle);
    }
};
class Cylinder : public MultiPolygonShape
{
  public:
    explicit Cylinder(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        multi_polygon_.addACircle(cylinder_center, cylinder_radius, 100, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class WeaklyCompressibleFluidInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit WeaklyCompressibleFluidInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
    };

    void update(size_t index_i, Real dt){};

  protected:
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> &rho_, &p_;
};

class FarFieldBoundary : public NonReflectiveBoundaryVariableCorrection
{
  public:
    explicit FarFieldBoundary(BaseInnerRelation &inner_relation) : NonReflectiveBoundaryVariableCorrection(inner_relation)
    {
        rho_farfield_ = rho0_f;
        sound_speed_ = c_f;
        vel_farfield_ = Vecd(U_f, 0.0);
    };
    virtual ~FarFieldBoundary(){};
};

#endif // EULERIAN_FLOW_AROUND_CYLINDER_H
