/**
 * @file 	fsi2.h
 * @brief 	This is the case file for the test of fluid - structure interaction.
 * @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#ifndef FSI2_CASE_H
#define FSI2_CASE_H

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 11.0;                         /**< Channel length. */
Real DH = 4.1;                          /**< Channel height. */
Real resolution_ref = 0.1;              /**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;         /**< Boundary width, determined by specific layer of boundary particles. */
Vec2d insert_circle_center(2.0, 2.0);   /**< Location of the cylinder center. */
Real insert_circle_radius = 0.5;        /**< Radius of the cylinder. */
/** Beam related parameters. */
Real bh = 0.4 * insert_circle_radius; /**< Height of the beam. */
Real bl = 7.0 * insert_circle_radius; /**< Length of the beam. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                            /**< Density. */
Real U_f = 1.0;                                               /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                        /**< Speed of sound. */
Real Re = 100.0;                                              /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 10.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e3;    /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

    return water_block_shape;
}
/** create a beam shape */
Real hbh = bh / 2.0;
Vec2d BLB(insert_circle_center[0], insert_circle_center[1] - hbh);
Vec2d BLT(insert_circle_center[0], insert_circle_center[1] + hbh);
Vec2d BRB(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] - hbh);
Vec2d BRT(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] + hbh);
// Beam observer location
StdVec<Vecd> beam_observation_location = {0.5 * (BRT + BRB)};
std::vector<Vecd> createBeamShape()
{
    std::vector<Vecd> beam_shape;
    beam_shape.push_back(BLB);
    beam_shape.push_back(BLT);
    beam_shape.push_back(BRT);
    beam_shape.push_back(BRB);
    beam_shape.push_back(BLB);

    return beam_shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
    outer_wall_shape.push_back(Vecd(-DL_sponge - BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, -BW));
    outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));

    return outer_wall_shape;
}
/**
 * @brief create inner wall shape
 */
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
    inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
    inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, DH));
    inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
    inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

    return inner_wall_shape;
}
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createBeamShape(), ShapeBooleanOps::sub);
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
class Insert : public MultiPolygonShape
{
  public:
    explicit Insert(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createBeamShape(), ShapeBooleanOps::add);
    }
};
/** create the beam base as constrain shape. */
MultiPolygon createBeamBaseShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(createBeamShape(), ShapeBooleanOps::sub);
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
class FluidObserverParticleGenerator : public ObserverParticleGenerator
{
  public:
    explicit FluidObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
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

#endif // FSI2_CASE_H
