/**
 * @file fsi2.h
 * @brief This is the case file for the test of fluid - structure interaction.
 * @details We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
 * @author Chi Zhang and Xiangyu Hu
 */

#ifndef FSI2_CASE_H
#define FSI2_CASE_H

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 50.0;                         /**< Channel length. */
Real DH = 5;                            /**< Channel height. */
Real resolution_ref = 0.3;              /**< Global reference resolution. */
Real DL_sponge = resolution_ref * 10.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;

//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.05;    /**< Density. */
Real U_f = 1.0;        /**< Characteristic velocity. */
Real c_f = 10.0 * U_f; /**< Speed of sound. */
Real Re = 100.0;       /**< Reynolds number. */
Real mu_f = rho0_f * U_f * 2 / Re;
Real Inlet_pressure = 2;
Real Outlet_pressure = 0.1;
//定义buffer的几何，虽然我不知道是干什么的
//========================================================================
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(DL - 2.5 * resolution_ref, 0.3 * DH);
/** inflow buffer parameters */
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 15);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
//=========================================================================
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 10.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e3;    /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//这两个东西后面好像不让改
Vec2d BRB(8, 3);
Vec2d BRT(10, 4);
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------

BoundingBox system_domain_bounds(Vec2d(0 - BW, -6 - BW), Vec2d(50 + BW, 15 + BW));
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a water block shape */
// 水的形状
std::vector<Vecd> water_block_shape{
    Vecd(0, 5), Vecd(20, 5), Vecd(25, 11), Vecd(35, 10),
    Vecd(50, 10), Vecd(50, 8), Vecd(45, 8), Vecd(33, 4), Vecd(25, 2), Vecd(40, -3),
    Vecd(50, -4), Vecd(50, -6), Vecd(40, -6), Vecd(20, -3), Vecd(0, 0), Vecd(0, 5)};

/** wall shape */
std::vector<Vecd> wall1_shape{
    Vecd(0, 5), Vecd(0, 5 + BW), Vecd(20, 5 + BW), Vecd(20, 5), Vecd(0, 5)};
std::vector<Vecd> wall2_shape{
    Vecd(0, 0), Vecd(20, -3), Vecd(20, -3 - BW), Vecd(0, 0 - BW), Vecd(0, 0)};
std::vector<Vecd> wall3_shape{
    Vecd(35, 10 + BW), Vecd(50, 10 + BW), Vecd(50, 10), Vecd(35, 10), Vecd(35, 10 + BW)};
std::vector<Vecd> wall4_shape{
    Vecd(45, 8), Vecd(50, 8), Vecd(50, 8 - BW), Vecd(45, 8 - BW), Vecd(45, 8)};
std::vector<Vecd> wall5_shape{
    Vecd(40, -3 + BW), Vecd(50, -4 + BW), Vecd(50, -4), Vecd(40, -3), Vecd(40, -3 + BW)};
std::vector<Vecd> wall6_shape{
    Vecd(40, -6), Vecd(50, -6), Vecd(50, -6 - BW), Vecd(40, -6 - BW), Vecd(40, -6)};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape{
    Vecd(0, 5 + BW), Vecd(20, 5 + BW), Vecd(25, 11 + BW), Vecd(35, 10 + BW),
    Vecd(50, 10 + BW), Vecd(50, 8 - BW), Vecd(45, 8 - BW), Vecd(33, 4 - BW), Vecd(27, 2), Vecd(40, -3 + BW),
    Vecd(50, -4 + BW), Vecd(50, -6 - BW), Vecd(40, -6 - BW), Vecd(20, -3 - BW), Vecd(0, 0 - BW), Vecd(0, 5 + BW)};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape{
    Vecd(0, 5), Vecd(20, 5), Vecd(25, 11), Vecd(35, 10),
    Vecd(50, 10), Vecd(50, 8), Vecd(45, 8), Vecd(33, 4), Vecd(25, 2), Vecd(40, -3),
    Vecd(50, -4), Vecd(50, -6), Vecd(40, -6), Vecd(20, -3), Vecd(0, 0), Vecd(0, 5)};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {

        multi_polygon_.addAPolygon(wall1_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(wall2_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(wall3_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(wall4_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(wall5_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(wall6_shape, ShapeBooleanOps::add);
    }
};

class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};
// 弹性部分
//弹性体是一整个管
class Insert : public MultiPolygonShape
{
  public:
    explicit Insert(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};//再把头尾的入口出口附近设置成约束不动
/** create the beam base as constrain shape. */
MultiPolygon createBeamBaseShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(wall1_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(wall2_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(wall3_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(wall4_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(wall5_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(wall6_shape, ShapeBooleanOps::add);

    return multi_polygon;
}
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        if (aligned_box_.checkInBounds(0, position))
        {
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        }
        return target_velocity;
    }
};
/** fluid observer particle generator */
class FluidObserver;
template <>
class ParticleGenerator<FluidObserver> : public ParticleGenerator<Observer>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body) : ParticleGenerator<Observer>(sph_body)
    {
        /** A line of measuring points at the entrance of the channel. */
        size_t number_observation_points = 21;
        Real range_of_measure = DH - resolution_ref * 4.0;
        Real start_of_measure = resolution_ref * 2.0;
        /** the measuring locations */
        for (size_t i = 0; i < number_observation_points; ++i)
        {
            Vec2d point_coordinate(0.0, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
            positions_.push_back(point_coordinate);
        }
    }
};
// namespace SPH
#endif // FSI2_CASE_H
