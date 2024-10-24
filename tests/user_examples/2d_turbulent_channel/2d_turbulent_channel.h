/**
 * @file 	2d_turbulent_channel.h
 * @brief 	Numerical parameters and body definition for 2d_turbulent_channel.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
//#include "k-epsilon_turbulent_model_complex.h"
//#include "k-epsilon_turbulent_model_complex.hpp"
#include "k-epsilon_turbulent_model.cpp"
using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Global parameters on the turbulent properties
//----------------------------------------------------------------------
Real y_p_theo = 0.05;      /**< Theoretical distance from the first particle P to wall  */
Real resolution_ref = 0.1; /**< Initial reference particle spacing. */
Real offset_dist_ref = y_p_theo - 0.5 * resolution_ref;
//Real offset_dist_ref = 0.0;
//** Initial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = {0.000180001, 3.326679e-5, 1.0e-9};
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 120;                        /**< Reference length. */
Real DH = 2 - 2.0 * offset_dist_ref;  /**< Reference and the height of main channel. */
Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
StdVec<int> id_exclude;
//-------------------------------------------------------
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + BW));
/** Observation locations, but for channel flow, the cell-based monitoring approach is used, parameters are defined in corresponding .cpp file*/
Real x_observe = 0.90 * DL;
Real x_observe_start = 0.90 * DL;
Real observe_spacing_x = 0.02 * DL;
int num_observer_points_x = 1;
int num_observer_points = std::round(DH / resolution_ref); //**Every particle is regarded as a cell monitor*
StdVec<Real> monitoring_bound = {109, 111};
Real observe_spacing = DH / num_observer_points;
StdVec<Vecd> observation_locations;
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;    /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f;
Real Re = 40000.0; /**< Reynolds number. */
//Real Re = 100.0;
Real mu_f = rho0_f * U_f * (DH + 2.0 * offset_dist_ref) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the emitter block. */
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d inlet_buffer_translation = Vec2d(-DL_sponge, 0.0) + inlet_buffer_halfsize;

/** the water block . */
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0),
    Vecd(-DL_sponge, DH),
    Vecd(DL, DH),
    Vecd(DL, 0.0),
    Vecd(-DL_sponge, 0.0)};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape{
    Vecd(-DL_sponge - 2.0 * BW, -BW),
    Vecd(-DL_sponge - 2.0 * BW, DH + BW),
    Vecd(DL + 2.0 * BW, DH + BW),
    Vecd(DL + 2.0 * BW, -BW),
    Vecd(-DL_sponge - 2.0 * BW, -BW),
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape{
    Vecd(-DL_sponge - 3.0 * BW, 0.0),
    Vecd(-DL_sponge - 3.0 * BW, DH),
    Vecd(DL + 3.0 * BW, DH),
    Vecd(DL + 3.0 * BW, 0.0),
    Vecd(-DL_sponge - 3.0 * BW, 0.0),
};
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
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
//----------------------------------------------------------------------
//	Define turbulent inflow boundary condition
//----------------------------------------------------------------------
