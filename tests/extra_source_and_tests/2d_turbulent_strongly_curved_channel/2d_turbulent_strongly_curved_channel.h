/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "k-epsilon_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.0254;
Real DH = 3.0 * scale; /**< Channel height. */
Real num_fluid_cross_section = 20.0;
Real central_angel = 90.0 * 2.0 * Pi / 360.0;
Real extend_in = 0.0;
Real extend_out = 0.0;
Real DL1 = 3.0 * DH + extend_in;
Real DL2 = 3.0 * DH + extend_out;
Real R1 = 27.0 * scale;
Real R2 = 30.0 * scale;

Vec2d circle_center(DL1, R2);
Real a = sin(central_angel);
Real b = cos(central_angel);
Vec2d point_N(R2 *a + DL1, R2 - R2 * b);
Vec2d point_M(R1 *a + DL1, R2 - R1 * b);
Vec2d vec_MN = point_N - point_M;
Vec2d vec_MN_vertical(-1.0 * vec_MN[1], vec_MN[0]);
Vec2d unit_vec_MN_vertical = vec_MN_vertical.normalized();
Vec2d point_P = point_N + DL2 * unit_vec_MN_vertical;
Vec2d point_Q = point_M + DL2 * unit_vec_MN_vertical;

Real DL_domain = point_P[0];
Real DH_domain = point_Q[1];
//----------------------------------------------------------------------
//	Unique parameters for turbulence.
//----------------------------------------------------------------------
Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
//** For K and Epsilon, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21 *
int type_turbulent_inlet = 0;
Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for AMRD *
int is_AMRD = 0;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
Real weight_vel_grad_sub_nearwall = 0.1;
//** Empirical parameter for initial stability*
Real turbulent_module_activate_time = 2.5;
//** Intial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = {0.000180001, 3.326679e-5, 1.0e-9};

Real y_p_constant = DH / 2.0 / num_fluid_cross_section; //** For the first try *
//Real y_p_constant = 0.05;
Real resolution_ref = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
Real offset_distance = y_p_constant - resolution_ref / 2.0;                        //** Basically offset distance is large than or equal to 0 *

Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;
Real half_channel_height = DH / 2.0;
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL_domain + 2.0 * BW, DH_domain + 2.0 * BW));

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
//** Laminar *
//Real U_inlet = 0.5;
//Real U_max = 0.75;
//Real U_f = U_inlet; //*Characteristic velo is regarded as average velo here
//Real c_f = 10.0 * U_max;                                        /**< Speed of sound. */
//Real rho0_f = 1.0;                                            /**< Density. */
//Real mu_f = 0.01;

//** Same parameter as SPH_4 *
Real U_inlet = 1.0;
Real Outlet_pressure = 0.0;
Real U_f = U_inlet;         //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */
Real Re = 148400.0;
//Real Re = 100.0;
Real mu_f = rho0_f * U_f * DH / Re;

Real Re_calculated = U_f * DH * rho0_f / mu_f;

Real DH_C = DH - 2.0 * offset_distance;

//** Intial inlet pressure to drive flow *
//Real initial_inlet_pressure = 0.5 * rho0_f * U_inlet ;
//Real initial_inlet_pressure = 100.0;
//----------------------------------------------------------------------
//	The emitter block with offset model.
//----------------------------------------------------------------------
Vec2d left_buffer_halfsize = Vec2d(0.5 * BW, 0.5 * DH_C + BW);
Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize + Vecd(0.0, offset_distance - BW);

Real outlet_buffer_length = BW;
Real outlet_buffer_height = 1.5 * DH;

//Real outlet_disposer_rotation_angel = 0.5 * Pi ; //** By default, counter-clockwise is positive *
//Vec2d outlet_buffer_center_translation = Vec2d(DL_domain - 0.5 * DH , DH_domain- 0.5 * outlet_buffer_length) ;

Real outlet_disposer_rotation_angel = central_angel; //** By default, counter-clockwise is positive *
Vec2d outlet_buffer_center_translation = (point_Q + point_P) / 2.0 - unit_vec_MN_vertical * outlet_buffer_length / 2.0;

Real outlet_emitter_rotation_angel = Pi + outlet_disposer_rotation_angel; //** By default, counter-clockwise is positive *

//** If return to the straight channel *
// Real outlet_disposer_rotation_angel = 0.0 * Pi ; //** By default, counter-clockwise is positive *
// Real outlet_emitter_rotation_angel =  Pi ; //** By default, counter-clockwise is positive *
// Vec2d outlet_buffer_center_translation = Vec2d(DL_domain - 0.5 * outlet_buffer_length , 0.5 * DH) ;

Vec2d right_buffer_halfsize = Vec2d(0.5 * outlet_buffer_length, 0.75 * outlet_buffer_height);
Vec2d right_buffer_translation = outlet_buffer_center_translation;

//Vec2d disposer_halfsize = Vec2d(0.75 * DH, 0.5 * BW);
//Vec2d disposer_translation = Vec2d(DL_domain + 0.25 * DH, DH_domain ) - disposer_halfsize;
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
// ** By kernel weight. *
int screen_output_interval = 100;
Real end_time = 10.0;              /**< End time. */
Real Output_Time = end_time / 2.0; /**< Time stamps for output of body states. */

Real cutoff_time = end_time * 0.6; //** cutoff_time should be a integral and the same as the PY script */
int number_observe_line = 4;
Real observe_angles[4] = {
    25.0 * (2.0 * Pi / 360.0),
    45.0 * (2.0 * Pi / 360.0),
    65.0 * (2.0 * Pi / 360.0),
    85.0 * (2.0 * Pi / 360.0)};
Real observer_offset_distance = 1.8 * resolution_ref;

int num_observer_points = std::round(DH_C / resolution_ref); //**Evrey particle is regarded as a cell monitor*
Real observe_spacing = DH_C / num_observer_points;

StdVec<Vecd> observation_locations;
void getPositionsOfMultipleObserveLines()
{
    for (int k = 0; k < number_observe_line; ++k)
    {
        Real ob_a = sin(observe_angles[k]);
        Real ob_b = cos(observe_angles[k]);
        Vec2d ob_point_at_outer_wall(R2 * ob_a + DL1, R2 - R2 * ob_b);
        Vec2d vec_point_at_outer_wall_to_circle_center = circle_center - ob_point_at_outer_wall;
        Vec2d unit_direction_observe = vec_point_at_outer_wall_to_circle_center.normalized();
        Vecd pos_observe_start = ob_point_at_outer_wall + y_p_constant * unit_direction_observe;
        for (int i = 0; i < num_observer_points; ++i)
        {
            Real offset = 0.0;
            offset = (i == 0 ? -observer_offset_distance : (i == num_observer_points - 1 ? observer_offset_distance : 0.0));
            Vecd pos_observer_i = pos_observe_start + (i * observe_spacing + offset) * unit_direction_observe;
            observation_locations.push_back(pos_observer_i);
        }
    }
}
void output_observe_positions()
{
    std::string filename = "../bin/output/observer_positions.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const Vecd &position : observation_locations)
    {
        outfile << position[0] << " " << position[1] << "\n";
    }
    outfile.close();
}
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
Real arc_sampling_interval = 0.01; //** Specify the arc resolution here *
int num_inner_arc_points = int(R1 * a / arc_sampling_interval);
int num_outer_arc_points = int(R2 * a / arc_sampling_interval);
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;

    //** 3 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, DH));
    water_block_shape.push_back(Vecd(DL1, DH));

    //** Inner Circle segment *
    Real start_x = DL1;
    for (int k = 1; k <= num_inner_arc_points; ++k)
    {
        Real x_coordinate = start_x + k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R1. Equation is (x-DL1)^2+(y-R2)^2=R1^2,
        //** Considring the coordinate, y= -sqr(R1^2-(x-DH1)^2)+R2 *
        Real y_coordinate = -sqrt(R1 * R1 - (x_coordinate - DL1) * (x_coordinate - DL1)) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** 4 points for outlet tube *
    water_block_shape.push_back(point_M);
    water_block_shape.push_back(point_Q + offset_distance * unit_vec_MN_vertical);
    water_block_shape.push_back(point_P + offset_distance * unit_vec_MN_vertical);
    water_block_shape.push_back(point_N);

    //** Outer Circle segment *
    start_x = point_N[0];
    for (int k = 1; k <= num_outer_arc_points; ++k)
    {
        Real x_coordinate = start_x - k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R2. Equation is (x-DL1)^2+(y-R2)^2=R2^2,
        //** Considring the coordinate, y= -sqr(R2^2-(x-DH1)^2)+R2 *
        Real y_coordinate = -sqrt(R2 * R2 - (x_coordinate - DL1) * (x_coordinate - DL1)) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** If return to straight channel, add extra 2 points *
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance, DH));
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance, 0.0));

    //** The left 2 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(DL1, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, 0.0));

    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance, computational_domain, "ComputationalDomain");
    }
};

std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> water_block_shape;

    //** 3 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - BW, DH));
    water_block_shape.push_back(Vecd(DL1, DH));

    //** Inner Circle segment *
    Real start_x = DL1;
    for (int k = 1; k <= num_inner_arc_points; ++k)
    {
        Real x_coordinate = start_x + k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R1. Equation is (x-DL1)^2+(y-R2)^2=R1^2,
        //** Considring the coordinate, y= -sqr(R1^2-(x-DH1)^2)+R2 *
        Real y_coordinate = -sqrt(R1 * R1 - (x_coordinate - DL1) * (x_coordinate - DL1)) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** 4 points for outlet tube *
    water_block_shape.push_back(point_M);
    water_block_shape.push_back(point_Q + (offset_distance + BW) * unit_vec_MN_vertical);
    water_block_shape.push_back(point_P + (offset_distance + BW) * unit_vec_MN_vertical);
    water_block_shape.push_back(point_N);

    //** Outer Circle segment *
    start_x = point_N[0];
    for (int k = 1; k <= num_outer_arc_points; ++k)
    {
        Real x_coordinate = start_x - k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R2. Equation is (x-DL1)^2+(y-R2)^2=R2^2,
        //** Considring the coordinate, y= -sqr(R2^2-(x-DH1)^2)+R2 *
        Real y_coordinate = -sqrt(R2 * R2 - (x_coordinate - DL1) * (x_coordinate - DL1)) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** If return to straight channel, add extra 2 points *
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + BW, DH));
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + BW, 0.0));

    //** The left 2 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(DL1, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - BW, 0.0));

    return water_block_shape;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> water_block_shape;

    //** 3 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - 2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - 2.0 * BW, DH));
    water_block_shape.push_back(Vecd(DL1, DH));

    //** Inner Circle segment *
    Real start_x = DL1;
    for (int k = 1; k <= num_inner_arc_points; ++k)
    {
        Real x_coordinate = start_x + k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R1. Equation is (x-DL1)^2+(y-R2)^2=R1^2,
        //** Considring the coordinate, y= -sqr(R1^2-(x-DH1)^2)+R2 *
        Real y_coordinate = -sqrt(R1 * R1 - (x_coordinate - DL1) * (x_coordinate - DL1)) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** 4 points for outlet tube *
    water_block_shape.push_back(point_M);
    water_block_shape.push_back(point_Q + (offset_distance + 2.0 * BW) * unit_vec_MN_vertical);
    water_block_shape.push_back(point_P + (offset_distance + 2.0 * BW) * unit_vec_MN_vertical);
    water_block_shape.push_back(point_N);

    //** Outer Circle segment *
    start_x = point_N[0];
    for (int k = 1; k <= num_outer_arc_points; ++k)
    {
        Real x_coordinate = start_x - k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R2. Equation is (x-DL1)^2+(y-R2)^2=R2^2,
        //** Considring the coordinate, y= -sqr(R2^2-(x-DH1)^2)+R2 *
        Real y_coordinate = -sqrt(R2 * R2 - (x_coordinate - DL1) * (x_coordinate - DL1)) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** If return to straight channel, add extra 2 points *
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + 2.0 * BW, DH));
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + 2.0 * BW, 0.0));

    //** The left 2 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(DL1, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - 2.0 * BW, 0.0));

    return water_block_shape;
}

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_dummy_boundary(createOuterWallShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance + BW, outer_dummy_boundary, "OuterDummyBoundary");

        MultiPolygon inner_dummy_boundary(createInnerWallShape());
        subtract<ExtrudeShape<MultiPolygonShape>>(-offset_distance, inner_dummy_boundary, "InnerDummyBoundary");
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
        : u_ref_(U_inlet), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        //target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        //target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_channel_height / half_channel_height);
        target_velocity[0] = u_ave;
        if (position[1] > half_channel_height)
        {
            std::cout << "Particles out of domain, wrong inlet velocity." << std::endl;
            std::cout << position[1] << std::endl;
            std::cin.get();
        }
        target_velocity[1] = 0.0;
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
        : Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_inlet), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(const Vecd &position) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
    }
};

struct RightOutflowPressure
{
    template <class BoundaryConditionType>
    RightOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        // Real run_time_ = GlobalStaticVariables::physical_time_;
        // if(run_time_ > 20.0)
        // {
        return p_;
        // }
        // else
        // {
        //     return initial_inlet_pressure;
        // }
    }
};