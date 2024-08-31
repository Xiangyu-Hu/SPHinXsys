/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "sphinxsys.h"
#include "k-epsilon_turbulent_model.cpp"
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 2.0;                         /**< Channel height. */
Real num_fluid_cross_section = 20.0;

Real extend_in = 0.0 ;
Real extend_out = 0.0 ;
Real DL1 = 4.0 + extend_in;
Real DL2 = 4.0 + extend_out;
Real R1 = 18.0;
Real R2 = R1 + DH;
Real DL_domain = DL1+R1+DH;
Real DH_domain = DH+R1+DL2;
Vec2d circle_center(DL1, R2);

//----------------------------------------------------------------------
//	Unique parameters for turbulence. 
//----------------------------------------------------------------------
Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
//** For K and Epsilon, type of the turbulent inlet, 0 is freestream, 1 is from interpolation from PY21 *
int type_turbulent_inlet = 0 ;
Real relaxation_rate_turbulent_inlet = 0.8;
//** Tag for AMRD *
int is_AMRD = 0 ;
//** Weight for correcting the velocity  gradient in the sub near wall region  *
Real weight_vel_grad_sub_nearwall_ = 0.1;
//** Intial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = { 0.000180001 ,3.326679e-5 ,1.0e-9 };  

Real y_p_constant = 0.05;
Real resolution_ref = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
Real offset_distance = y_p_constant - resolution_ref / 2.0; //** Basically offset distance is large than or equal to 0 *

Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
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
Real U_f = U_inlet; //*Characteristic velocity
Real U_max = 1.5 * U_inlet; //** An estimated value, generally 1.5 U_inlet *
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0;                                            /**< Density. */
Real Re = 40000.0;
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
Vec2d left_buffer_halfsize = Vec2d(0.5 * BW, 0.5 * DH_C + BW );
Vec2d left_buffer_translation = Vec2d(-DL_sponge, 0.0) + left_buffer_halfsize + Vecd(0.0, offset_distance - BW);

Real outlet_buffer_length = BW ;
Real outlet_buffer_height = 1.5 * DH ;

Real outlet_disposer_rotation_angel = 0.5 * Pi ; //** By default, counter-clockwise is positive *
Real outlet_emitter_rotation_angel = -0.5 * Pi ; //** By default, counter-clockwise is positive *
Vec2d outlet_buffer_center_translation = Vec2d(DL_domain - 0.5 * DH , DH_domain- 0.5 * outlet_buffer_length) ;

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
Real x_observe_start = 0.99 * DL1;
int num_observer_points = std::round(DH_C / resolution_ref); //**Evrey particle is regarded as a cell monitor* 
Real observe_spacing = DH_C / num_observer_points;

// By kernel weight.
StdVec<Vecd> observation_location;
Vecd pos_observe_start = Vecd(x_observe_start, resolution_ref/2.0 + offset_distance);
Vecd unit_direction_observe = Vecd(0.0, 1.0);
Real observer_offset_distance = 2.0 * resolution_ref ;
//----------------------------------------------------------------------
//	Cases-dependent geometries 
//----------------------------------------------------------------------
Real arc_sampling_interval = 0.01;   //** Specify the arc resolution here *
int num_inner_arc_points = int(R1 / arc_sampling_interval);
int num_outer_arc_points = int(R2 / arc_sampling_interval);
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    
    //** 3 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, DH));
    water_block_shape.push_back(Vecd(DL1, DH));
    

    // //** Inner Circle segment *
    // Real start_x = DL1;
    // for (int k = 1; k <= num_inner_arc_points; ++k)
    // {
    //     Real x_coordinate = start_x + k * arc_sampling_interval; //** clockwise *
    //     //** Circle center is (DL1, R2), radius is R1. Equation is (x-DL1)^2+(y-R2)^2=R1^2,
    //     //** Considring the coordinate, y= -sqr(R1^2-(x-DH1)^2)+R2 * 
    //     Real y_coordinate = -sqrt( R1*R1 - (x_coordinate-DL1)*(x_coordinate-DL1) ) + R2;
    //     water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    // }

    // //** 4 points for outlet tube *
    // water_block_shape.push_back(Vecd(DL1 + R1, R2));
    // water_block_shape.push_back(Vecd(DL1 + R1, R2 + DL2 + offset_distance));
    // water_block_shape.push_back(Vecd(DL1 + R2, R2 + DL2 + offset_distance));
    // water_block_shape.push_back(Vecd(DL1 + R2, R2 ));

    // //** Outer Circle segment *
    // start_x = DL1 + R2;
    // for (int k = 1; k <= num_outer_arc_points; ++k)
    // {
    //     Real x_coordinate = start_x - k * arc_sampling_interval; //** clockwise *
    //     //** Circle center is (DL1, R2), radius is R2. Equation is (x-DL1)^2+(y-R2)^2=R2^2,
    //     //** Considring the coordinate, y= -sqr(R2^2-(x-DH1)^2)+R2 * 
    //     Real y_coordinate = -sqrt( R2*R2 - (x_coordinate-DL1)*(x_coordinate-DL1) ) + R2;
    //     water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    // }

    //** If return to straight channel, add extra 2 points *
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance, DH));
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance, 0.0));

    //** The left 2 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(DL1, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance,  0.0));
    
    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
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
        Real y_coordinate = -sqrt( R1*R1 - (x_coordinate-DL1)*(x_coordinate-DL1) ) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** 4 points for outlet tube *
    water_block_shape.push_back(Vecd(DL1 + R1, R2));
    water_block_shape.push_back(Vecd(DL1 + R1, R2 + DL2 + offset_distance + BW));
    water_block_shape.push_back(Vecd(DL1 + R2, R2 + DL2 + offset_distance + BW));
    water_block_shape.push_back(Vecd(DL1 + R2, R2 ));

    //** Outer Circle segment *
    start_x = DL1 + R2;
    for (int k = 1; k <= num_outer_arc_points; ++k)
    {
        Real x_coordinate = start_x - k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R2. Equation is (x-DL1)^2+(y-R2)^2=R2^2,
        //** Considring the coordinate, y= -sqr(R2^2-(x-DH1)^2)+R2 * 
        Real y_coordinate = -sqrt( R2*R2 - (x_coordinate-DL1)*(x_coordinate-DL1) ) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** If return to straight channel, add extra 2 points *
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + BW, DH));
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + BW, 0.0));

    //** The left 2 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(DL1, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - BW,  0.0));
    
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
        Real y_coordinate = -sqrt( R1*R1 - (x_coordinate-DL1)*(x_coordinate-DL1) ) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** 4 points for outlet tube *
    water_block_shape.push_back(Vecd(DL1 + R1, R2));
    water_block_shape.push_back(Vecd(DL1 + R1, R2 + DL2 + offset_distance + 2.0 * BW));
    water_block_shape.push_back(Vecd(DL1 + R2, R2 + DL2 + offset_distance + 2.0 * BW));
    water_block_shape.push_back(Vecd(DL1 + R2, R2 ));

    //** Outer Circle segment *
    start_x = DL1 + R2;
    for (int k = 1; k <= num_outer_arc_points; ++k)
    {
        Real x_coordinate = start_x - k * arc_sampling_interval; //** clockwise *
        //** Circle center is (DL1, R2), radius is R2. Equation is (x-DL1)^2+(y-R2)^2=R2^2,
        //** Considring the coordinate, y= -sqr(R2^2-(x-DH1)^2)+R2 * 
        Real y_coordinate = -sqrt( R2*R2 - (x_coordinate-DL1)*(x_coordinate-DL1) ) + R2;
        water_block_shape.push_back(Vecd(x_coordinate, y_coordinate));
    }

    //** If return to straight channel, add extra 2 points *
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + 2.0 * BW, DH));
    // water_block_shape.push_back(Vecd(DL_domain + offset_distance + 2.0 * BW, 0.0));


    //** The left 2 points for inlet tube, totally 5 *
    water_block_shape.push_back(Vecd(DL1, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance - 2.0 * BW,  0.0));
    
    return water_block_shape;
}

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_dummy_boundary(createOuterWallShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance + BW, outer_dummy_boundary, "OuterDummyBoundary");

        MultiPolygon inner_dummy_boundary(createInnerWallShape());
        subtract<ExtrudeShape<MultiPolygonShape>>(-offset_distance, inner_dummy_boundary, "InnerDummyBoundary");
    }
};
// class WallBoundary : public MultiPolygonShape
// {
// public:
//     explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
//     {
//         /** Geometry definition. */
//         //** Get the ring *
//         multi_polygon_.addACircle(circle_center, R2+BW, 100, ShapeBooleanOps::add);
//         multi_polygon_.addACircle(circle_center, R2, 100, ShapeBooleanOps::sub);
//         multi_polygon_.addACircle(circle_center, R1 , 100, ShapeBooleanOps::add);
//         multi_polygon_.addACircle(circle_center, R1-BW, 100, ShapeBooleanOps::sub);
//         //** Get the 1/4 ring *
//         std::vector<Vecd> sub_circle;
//         Vec2d temp =  Vec2d::Zero();
//         temp = circle_center + Vec2d(-(R2 + BW), -(R2 + BW));
//         sub_circle.push_back(temp);
//         temp = circle_center + Vec2d(-(R2 + BW), (R2 + BW));
//         sub_circle.push_back(temp);
//         temp = circle_center + Vec2d((R2 + BW), (R2 + BW));
//         sub_circle.push_back(temp);
//         temp = circle_center + Vec2d((R2 + BW), 0.0);
//         sub_circle.push_back(temp);
//         temp = circle_center + Vec2d(0.0, 0.0);
//         sub_circle.push_back(temp);
//         temp = circle_center + Vec2d(0.0, -(R2 + BW));
//         sub_circle.push_back(temp);
//         temp = circle_center + Vec2d(-(R2 + BW), -(R2 + BW));
//         sub_circle.push_back(temp);
//         multi_polygon_.addAPolygon(sub_circle, ShapeBooleanOps::sub);

//         //** Get the inlet straight channel wall *
//         std::vector<Vecd> outer_wall_shape1;
//         outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -BW));
//         outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, DH + BW));
//         outer_wall_shape1.push_back(Vecd(DL1 , DH + BW));
//         outer_wall_shape1.push_back(Vecd(DL1 , -BW));
//         outer_wall_shape1.push_back(Vecd(-DL_sponge - BW, -BW));
//         multi_polygon_.addAPolygon(outer_wall_shape1, ShapeBooleanOps::add);

//         std::vector<Vecd> inner_wall_shape1;
//         inner_wall_shape1.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0 ));
//         inner_wall_shape1.push_back(Vecd(-DL_sponge - 2.0 * BW, DH ));
//         inner_wall_shape1.push_back(Vecd(DL1, DH ));
//         inner_wall_shape1.push_back(Vecd(DL1, 0.0));
//         inner_wall_shape1.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0 ));
//         multi_polygon_.addAPolygon(inner_wall_shape1, ShapeBooleanOps::sub);
        
//         //** Get the outlet straight channel wall *
//         std::vector<Vecd> outer_wall_shape2;
//         outer_wall_shape2.push_back(Vecd(DL1 + R1 - BW     , DH + R1));
//         outer_wall_shape2.push_back(Vecd(DL1 + R1 - BW     , DH + R1 + DL2 + BW ));
//         outer_wall_shape2.push_back(Vecd(DL1 + R1+ DH + BW , DH + R1 + DL2 + BW ));
//         outer_wall_shape2.push_back(Vecd(DL1 + R1+ DH + BW , DH + R1 ));
//         outer_wall_shape2.push_back(Vecd(DL1 + R1 - BW     , DH + R1));
//         multi_polygon_.addAPolygon(outer_wall_shape2, ShapeBooleanOps::add);

//         std::vector<Vecd> inner_wall_shape2;
//         inner_wall_shape2.push_back(Vecd(DL1 + R1      , DH + R1));
//         inner_wall_shape2.push_back(Vecd(DL1 + R1      , DH + R1 + DL2+ 2.0 * BW));
//         inner_wall_shape2.push_back(Vecd(DL1 + R1 + DH , DH + R1 + DL2 + +2.0 * BW));
//         inner_wall_shape2.push_back(Vecd(DL1 + R1 + DH , DH + R1));
//         inner_wall_shape2.push_back(Vecd(DL1 + R1      , DH + R1));
//         multi_polygon_.addAPolygon(inner_wall_shape2, ShapeBooleanOps::sub);
//     }
// };
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape& aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType& boundary_condition)
        : u_ref_(U_inlet), t_ref_(2.0),
        aligned_box_(boundary_condition.getAlignedBox()),
        halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd& position, Vecd& velocity)
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
            std::cout <<position[1]<<std::endl;
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
    Real  t_ref_, u_ref_, du_ave_dt_;

public:
    explicit TimeDependentAcceleration(Vecd gravity_vector)
        : Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_inlet), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(const Vecd& position) override
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