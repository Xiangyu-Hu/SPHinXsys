/**
 * @file 	2d_laminar_wavy_channel.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "sphinxsys.h"
#include "k-epsilon_turbulent_model_complex.h"
#include "k-epsilon_turbulent_model_complex.hpp"
#include "k-epsilon_turbulent_model.cpp"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 8.0;                         /**< Channel length. */
Real DH = 1.0;                         /**< Channel height. */
Real resolution_ref = 0.05;              /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 10.0; /**< Sponge region to impose inflow condition. */
Real amplitude = 0.1;
Real wave_length = 1;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge-3.0*BW, -3.0 * BW), Vec2d(DL + 3.0 * BW, DH + 3.0 * BW));
// Observation locations
Vec2d point_coordinate_1(3.0, 5.0);
Vec2d point_coordinate_2(4.0, 5.0);
Vec2d point_coordinate_3(5.0, 5.0);
StdVec<Vec2d> observation_locations = {point_coordinate_1, point_coordinate_2, point_coordinate_3};
//----------------------------------------------------------------------
//	Global parameters on the turbulent properties
//----------------------------------------------------------------------
//Real y_p_theo = 0.05;                 /**< Theoretical distance from the first particle P to wall  */
Real y_p_theo = 0.0;
//Real offset_dist_ref = y_p_theo - 0.5 * resolution_ref;
Real offset_dist_ref = 0.0;
StdVec<int> id_exclude;
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                            /**< Density. */
Real U_f = 0.816;                                               /**< freestream velocity. */
//Real U_f = 1.0;
Real c_f = 10.0 * U_f;                                        /**< Speed of sound. */
//Real Re = 100.0;                                              /**< Reynolds number. */
//Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//Real mu_f = 0.0001; /**< Dynamics viscosity. */
Real mu_f = 0.01; /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH ));
    water_block_shape.push_back(Vecd(DL, DH ));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

    return water_block_shape;
}
/** create a water block buffer shape. */
MultiPolygon createBufferShape()
{
    std::vector<Vecd> buffer_shape;
    buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
    buffer_shape.push_back(Vecd(-DL_sponge, DH ));
    buffer_shape.push_back(Vecd(0.0, DH));
    buffer_shape.push_back(Vecd(0.0, 0.0));
    buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(buffer_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
/** the water in the trough part. */
std::vector<Vecd> water_shape_trough
{
    Vecd(0.0, 0.0),
    Vecd(DL, 0.0),
    Vecd(DL, -amplitude),
    Vecd(0.0, -amplitude),
    Vecd(0.0, 0.0),
};
/** the top wall polygon. */
std::vector<Vecd> top_wall_shape
{
    Vecd(-DL_sponge - BW, DH),
    Vecd(-DL_sponge - BW, DH + BW),
    Vecd(DL + BW, DH + BW),
    Vecd(DL + BW, DH),
    Vecd(-DL_sponge - BW, DH),
};
/** the bottom wall polygon. */
std::vector<Vecd> bottom_wall_shape{};
class WavyShape
{
    int num_wavy_points;
    Real wavy_interval = 0.01;
public:
    explicit WavyShape()
    {
        bottom_wall_shape.push_back(Vecd(-DL_sponge - BW, 0.0));
        bottom_wall_shape.push_back(Vecd(0.0, 0.0));
        num_wavy_points = (DL - 0.0) / wavy_interval;
        for (size_t k = 1; k != num_wavy_points; ++k)
        {
            Real x_coordinate = k * wavy_interval;
            Real y_coordinate = -1.0 * amplitude * std::sin(2.0 * M_PI * x_coordinate);
            bottom_wall_shape.push_back(Vecd(x_coordinate,y_coordinate));
        }
        bottom_wall_shape.push_back(Vecd(DL, 0.0));
        bottom_wall_shape.push_back(Vecd(DL + BW, 0.0));
        bottom_wall_shape.push_back(Vecd(DL + BW, -BW - amplitude));
        bottom_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW - amplitude));
        bottom_wall_shape.push_back(Vecd(-DL_sponge - BW, 0.0));
    }
};
class ChannelShape
{
public:
    explicit ChannelShape()
    {
        bottom_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
        bottom_wall_shape.push_back(Vecd(-DL_sponge - BW, 0.0));
        bottom_wall_shape.push_back(Vecd(DL + BW, 0.0));
        bottom_wall_shape.push_back(Vecd(DL + BW, -BW));
        bottom_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
    }
};
/** Water block shape definition */
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        MultiPolygon multi_polygon;
        multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon.addAPolygon(water_shape_trough, ShapeBooleanOps::add);
        multi_polygon.addAPolygon(bottom_wall_shape, ShapeBooleanOps::sub);
        add<MultiPolygonShape>(multi_polygon);
    }
};

//----------------------------------------------------------------------
//	Definition of the wall
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(top_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(bottom_wall_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	Case dependent flow boundary condition.
//----------------------------------------------------------------------
class FreeStreamCondition : public fluid_dynamics::FlowVelocityBuffer
{
    Real u_ave_, u_ref_, t_ref;

  public:
    FreeStreamCondition(BodyPartByCell &constrained_region)
        : fluid_dynamics::FlowVelocityBuffer(constrained_region),
          u_ave_(0), u_ref_(U_f), t_ref(2.0) {}
    Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
    {
        /* Fully-developed velocity inlet */
        //Real radius = 0.0;
        //if (position[1]>= DH/2.0 )
        //{
        //    radius = position[1] - DH / 2.0;
        //}
        //else if (position[1] < DH / 2.0)
        //{
        //    radius = DH / 2.0 - position[1] ;
        //}
        //u_ave_ = 1.5 * u_ref_ * (1.0 - radius * radius / (DH/2.0) / (DH / 2.0));
        return Vecd(u_ave_, 0.0);
    }
    void setupDynamics(Real dt = 0.0) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;
    }
};
class CorrectBufferVelocity : public fluid_dynamics::BaseFlowBoundaryCondition
{
public:
    CorrectBufferVelocity(BodyPartByCell& body_part) : BaseFlowBoundaryCondition(body_part) {}
    virtual ~CorrectBufferVelocity() {};
    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i][1] = 0.0;
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
        : Gravity(gravity_vector), t_ref_(2.0), u_ref_(2.0*U_f), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(const Vecd& position) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
    }
};

