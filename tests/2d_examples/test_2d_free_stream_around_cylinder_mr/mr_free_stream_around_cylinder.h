/**
 * @file	mr_freestream_flow_around_cylinder.h
 * @brief 	We consider a flow pass the cylinder with freestream boundary condition in 2D.
 * @details Particles with different sizes are used in this test.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 30.0;                               /**< Channel length. */
Real DH = 16.0;                               /**< Channel height. */
Real particle_spacing_ref = 0.4;              /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;         /**< Sponge region to impose injection. */
Vec2d insert_circle_center(10.0, 0.5 * DH);   /**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;              /**< Radius of the cylinder. */
// Observation locations
Vec2d point_coordinate_1(3.0, 5.0);
Vec2d point_coordinate_2(4.0, 5.0);
Vec2d point_coordinate_3(5.0, 5.0);
StdVec<Vecd> observation_locations = {point_coordinate_1, point_coordinate_2, point_coordinate_3};
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                            /**< Density. */
Real U_f = 1.0;                                               /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                        /**< Speed of sound. */
Real Re = 100.0;                                              /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(-DL_sponge, 0.0)};
std::vector<Vecd> initial_refinement_region{
    Vecd(-DL_sponge - BW, 5.0), Vecd(-DL_sponge - BW, 11.0), Vecd(DL + BW, 11.0), Vecd(DL + BW, 5.0), Vecd(-DL_sponge - BW, 5.0)};
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;

//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(water_block_shape);
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon circle(insert_circle_center, insert_circle_radius, 100);
        subtract<MultiPolygonShape>(circle);
    }
};
class Cylinder : public MultiPolygonShape
{
  public:
    explicit Cylinder(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
struct FreeStreamVelocity
{
    Real u_ref_, t_ref_;

    template <class BoundaryConditionType>
    FreeStreamVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;
        target_velocity[0] = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
    Real t_ref_, u_ref_, du_ave_dt_;

  public:
    explicit TimeDependentAcceleration(Vecd gravity_vector)
        : Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_f), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(Vecd &position) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);

        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
    }
};
