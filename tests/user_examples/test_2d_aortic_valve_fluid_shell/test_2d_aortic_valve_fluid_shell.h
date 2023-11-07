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
Real scale = 0.001;
Real DH = 20 * scale;                   /**< Channel height. */
Real DL = 5 * DH;                       /**< Channel length. */
Real circle_radius = DH;                /**< Radius of the circle. */
Vecd circle_center(2.5 * DH, DH);       /**< Position of the circle center. */
Real resolution_ref = 0.5 * scale;      /**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;         /**< Boundary width, determined by specific layer of boundary particles. */
/** Beam related parameters. */
Real beam_position = circle_center[0] - circle_radius; /**< Position of beam fixed end. */
Real beam_angle = 45;                                  /**< Angle between beam and verticle position. */
Real thickness = 0.3 * scale;                          /**< Thickeness of the beam. */
Real bl = 26 * scale;                                  /**< Length of the beam. */
Real width = 117 * scale;                              /**<Width out of the plane>*/
Real resolution_beam = thickness;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + circle_radius + BW));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1000.0;           /**< Density. */
Real mu_f = 4.3e-3;             /**< Dynamics viscosity. */
Real Q_max = 23.44 * 1e-3 / 60; /**< Maximum flow rate in a period. */
Real U_f = Q_max / DH / width;  /**< Characteristic velocity. */
Real U_max = 3 * U_f;
Real c_f = 10.0 * U_max;            /**< Speed of sound. */
Real Re = rho0_f * U_f * DH / mu_f; /**< Reynolds number. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1000.0;        /**< Reference density.*/
Real poisson = 0.49;         /**< Poisson ratio.*/
Real Youngs_modulus = 1.5e6; /**< Youngs Modulus. */
Real shape_constant = 0.4;
Real physical_viscosity = shape_constant / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * thickness;
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
Real hbh = thickness / 2.0;
Real cos_a = cos(beam_angle);
Real sin_a = sin(beam_angle);
Vec2d BLB(beam_position + bl * sin_a - hbh * cos_a, DH - bl * cos_a - hbh * sin_a);
Vec2d BLT(beam_position - hbh * cos_a, DH - hbh * sin_a);
Vec2d BRT(beam_position + hbh * cos_a, DH + hbh * sin_a);
Vec2d BRB(beam_position + bl * sin_a + hbh * cos_a, DH - bl * cos_a + hbh * sin_a);
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
        multi_polygon_.addACircle(circle_center, circle_radius, 100, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createBeamShape(), ShapeBooleanOps::sub);
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(circle_center, circle_radius + BW, 100, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addACircle(circle_center, circle_radius, 100, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
/** Particle generator and constraint boundary for shell baffle. */
int particle_number_mid_surface = int(bl / resolution_beam);
class BeamParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit BeamParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            Real x = beam_position + (Real(i) + 0.5) * resolution_beam * sin_a;
            Real y = DH - (Real(i) + 0.5) * resolution_beam * cos_a;
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_beam);
            Vec2d normal_direction = Vec2d(cos_a, sin_a);
            initializeSurfaceProperties(normal_direction, thickness);
        }
    }
};
/** create the beam base as constrain shape. */
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry(){};

  private:
    void tagManually(size_t index_i)
    {
        Real dx = base_particles_.pos_[index_i][0] - beam_position;
        Real dy = base_particles_.pos_[index_i][1] - DH;
        Real distance = dx * dx + dy * dy;
        if (distance <= 16 * resolution_beam * resolution_beam)
            body_part_particles_.push_back(index_i);
    };
};
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
        : u_ref_(U_f), t_ref_(1.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / (t_ref_ / 2)));
        target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
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
