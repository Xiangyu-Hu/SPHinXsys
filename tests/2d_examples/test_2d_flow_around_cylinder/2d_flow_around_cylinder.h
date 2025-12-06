/**
 * @file 	2d_cylinder_flow.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 15.0;                         /**< Channel length. */
Real DH = 10.0;                         /**< Channel height. */
Real resolution_ref = 0.2;              /**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 10.0; /**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0;  /**< Sponge region to impose freestream condition. */
Vec2d insert_circle_center(4.0, 5.0);   /**< Location of the cylinder center. */
Real insert_circle_radius = 0.75;       /**< Radius of the cylinder. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
// Observation locations
Vec2d point_coordinate_1(3.0, 5.0);
Vec2d point_coordinate_2(4.0, 5.0);
Vec2d point_coordinate_3(5.0, 5.0);
StdVec<Vec2d> observation_locations = {point_coordinate_1, point_coordinate_2, point_coordinate_3};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                            /**< Density. */
Real U_f = 1.0;                                               /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                        /**< Speed of sound. */
Real Re = 100.0;                                              /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

    return water_block_shape;
}
/** create a water block buffer shape. */
MultiPolygon createBufferShape()
{
    std::vector<Vecd> buffer_shape;
    buffer_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
    buffer_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
    buffer_shape.push_back(Vecd(DL, DH + DH_sponge));
    buffer_shape.push_back(Vecd(DL, DH));
    buffer_shape.push_back(Vecd(0.0, DH));
    buffer_shape.push_back(Vecd(0.0, 0.0));
    buffer_shape.push_back(Vecd(DL, 0.0));
    buffer_shape.push_back(Vecd(DL, -DH_sponge));
    buffer_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(buffer_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
/** Water block shape definition */
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        MultiPolygon multi_polygon;
        multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
        add<MultiPolygonShape>(multi_polygon);
    }
};
//----------------------------------------------------------------------
//	Define parametrization for this case.
//----------------------------------------------------------------------
class ParameterizedViscosity : public BaseParameterization<Viscosity>
{
  public:
    ParameterizedViscosity(ConstructArgs<ParameterizationIO *, Real> args)
        : BaseParameterization<Viscosity>(std::get<0>(args), std::get<1>(args))
    {
        getAParameter("WaterMaterial", "Viscosity", mu_);
    }
};
//----------------------------------------------------------------------
//	Definition of the cylinder
//----------------------------------------------------------------------
class Cylinder : public MultiPolygonShape
{
  public:
    explicit Cylinder(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Case dependent flow boundary condition.
//----------------------------------------------------------------------
class FreeStreamCondition : public fluid_dynamics::FlowVelocityBuffer
{
    Real u_ave_, u_ref_, t_ref;
    Real *physical_time_;

  public:
    FreeStreamCondition(BodyPartByCell &constrained_region)
        : fluid_dynamics::FlowVelocityBuffer(constrained_region),
          u_ave_(0), u_ref_(U_f), t_ref(2.0),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {}
    Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
    {
        return Vecd(u_ave_, 0.0);
    }
    void setupDynamics(Real dt = 0.0) override
    {
        Real run_time = *physical_time_;
        u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;
    }
};
