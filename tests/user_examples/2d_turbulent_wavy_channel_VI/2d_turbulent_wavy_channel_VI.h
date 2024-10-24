/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "k-epsilon_turbulent_model.cpp"
#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"
#include "sphinxsys.h"
using namespace SPH;

/**
 * @ Parameters for validation with [1997 Joseph JCP].
 */
/*Real DH = 1.0e-3;
Real DL = 4.0 * DH;
Real rho0_f = 1000.0;
Real gravity_g = 1.0e-4;
Real nu_f = 1.0e-6;
Real mu_f = nu_f * rho0_f;
Real U_f = 1.25e-5;
Real c_f = 10.0 * U_f;*/

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 1.0;  /**< Channel height. */
Real DL = 40.0; /**< Channel length. */
Real amplitude = 0.1;
Real wave_length = 1;
Real resolution_ref = DH / 20.0; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;    /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;
Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
Real extend_out = 0.1 * DL;
Real extend_in = 0.1 * DL;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -DH), Vec2d(DL + 2.0 * BW + extend_out + extend_in, DH + 2.0 * BW));
StdVec<int> id_exclude;
Real y_p_theo = 0.0;
Real offset_dist_ref = 0.0;
//** Initial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = {0.01, 0.1, 0.001};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_f = 0.816;      //*Characteristic velo is regarded as average velo here
Real c_f = 10.0 * U_f; /**< Speed of sound. */
Real rho0_f = 1.0;     /**< Density. */
Real mu_f = 0.0001;
Real Re = U_f * DH * rho0_f / mu_f;
//----------------------------------------------------------------------
// Observation.
//----------------------------------------------------------------------
//** Trough *
int num_observer_points = std::round((DH + amplitude) / resolution_ref); //**Every particle is regarded as a cell monitor*
Real pos_x_cross_line = 1.25;
Real size_cell_x = 1.8 * resolution_ref;
StdVec<Real> monitoring_bound = {pos_x_cross_line - size_cell_x / 2.0, pos_x_cross_line + size_cell_x / 2.0}; //**The 2nd periodic wave*
Real observe_offset_y = -1.0 * amplitude;
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
/** the emitter block. */
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d inlet_buffer_translation = Vec2d(-DL_sponge, 0.0) + inlet_buffer_halfsize;
Vec2d disposer_up_halfsize = Vec2d(0.5 * BW, 0.55 * DH);
Vec2d disposer_up_translation = Vec2d(DL + extend_out + extend_in, -0.05 * DH) + disposer_up_halfsize;

/** the water block . */
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH),
    Vecd(DL + extend_out + extend_in, DH), Vecd(DL + extend_out + extend_in, 0.0), Vecd(-DL_sponge, 0.0)};
/** the water in the trough part. */
std::vector<Vecd> water_shape_trough{
    Vecd(0.0 + extend_in, 0.0),
    Vecd(DL + extend_in, 0.0),
    Vecd(DL + extend_in, -amplitude),
    Vecd(0.0 + extend_in, -amplitude),
    Vecd(0.0 + extend_in, 0.0),
};
/** the top wall polygon. */
std::vector<Vecd> top_wall_shape{
    Vecd(-DL_sponge - 2.0 * BW, DH),
    Vecd(-DL_sponge - 2.0 * BW, DH + BW),
    Vecd(DL + 2.0 * BW + extend_out + extend_in, DH + BW),
    Vecd(DL + 2.0 * BW + extend_out + extend_in, DH),
    Vecd(-DL_sponge - 2.0 * BW, DH),
};
/** the bottom wall polygon. */
std::vector<Vecd> bottom_wall_shape{};
class WavyShape
{
    int num_wavy_points, num_wavy_extend_left, num_wavy_extend_right;
    int num_total;
    Real wavy_interval = 0.01; //** Specify the wavy resolution here *
  public:
    explicit WavyShape()
    {
        bottom_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
        bottom_wall_shape.push_back(Vecd(0.0, 0.0));
        num_wavy_points = (DL - 0.0) / wavy_interval;
        for (int k = 1; k <= num_wavy_points; ++k)
        {
            Real x_coordinate = k * wavy_interval + extend_in;
            Real y_coordinate = -1.0 * amplitude * std::sin(2.0 * M_PI * x_coordinate);
            bottom_wall_shape.push_back(Vecd(x_coordinate, y_coordinate));
        }
        bottom_wall_shape.push_back(Vecd(DL + extend_in, 0.0));
        bottom_wall_shape.push_back(Vecd(DL + 2.0 * BW + extend_out + extend_in, 0.0));
        bottom_wall_shape.push_back(Vecd(DL + 2.0 * BW + extend_out + extend_in, -BW - amplitude));
        bottom_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, -BW - amplitude));
        bottom_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
    }
};

//----------------------------------------------------------------------
// @brief 	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(water_shape_trough, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(bottom_wall_shape, ShapeBooleanOps::sub);
    }
};

//----------------------------------------------------------------------
// @brief 	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        multi_polygon_.addAPolygon(top_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(bottom_wall_shape, ShapeBooleanOps::add);
    }
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
            /* Fully-developed velocity inlet */
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);

            /* Uniform velocity inlet */
            //target_velocity[0] = u_ave;

            /* Fix velocity in Y direction */
            target_velocity[1] = 0.0;
        }
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
    Real du_ave_dt_, u_ref_, t_ref_;

  public:
    explicit TimeDependentAcceleration(Vecd gravity_vector)
        : Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_f), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(const Vecd &position) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
    }
};