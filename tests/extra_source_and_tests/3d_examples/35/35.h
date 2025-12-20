/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
//#include "k-epsilon_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
#include "zeroth-order_residue.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
//% Dimension: 10 \mu m, s, kg
Real scale = 1.0e-6;
Real H_inlet = 4.0;
Real L_inlet = 10.0;
Real Radius_chamber = 83.0 / 2.0;
Real H_outlet = H_inlet;
Real L_outlet = L_inlet;
Real H_total = 100.0;

Real D_hydraulic = 2.0 * L_inlet * H_inlet / (L_inlet + H_inlet);

Real num_fluid_cross_H_inlet = 5.0;
Real resolution_ref = H_inlet / num_fluid_cross_H_inlet; /**< Initial reference particle spacing. */

Real BW = resolution_ref * 4; /**< Reference size of the emitter. */
Real buffer_thickness = 5.0 * resolution_ref;

//** STL relevant parameters */
Vecd point_O(0.0, 0.0, 0.0);
Vecd point_A = point_O + Vecd(0.0, 0.0, H_total);

Vecd point_OA_half = (point_O + point_A) / 2.0;

Real length_outlet = L_inlet; //% Length of the outlet channel, Geometry dependent
//----------------------------------------------------------------------
//	Domain bounds of the system, STL relevant.
//----------------------------------------------------------------------
Real extend_domain_length = length_outlet; //% Geometry dependent
BoundingBoxd system_domain_bounds(point_O +
                                     Vecd(-Radius_chamber, -Radius_chamber, 0.0) +
                                     Vecd(-length_outlet, -length_outlet, 0.0) +
                                     2.0 * Vecd(-BW, -BW, -BW) +
                                     Vecd(-extend_domain_length, -extend_domain_length, 0.0),
                                 point_A +
                                     Vecd(Radius_chamber, Radius_chamber, 0.0) +
                                     Vecd(length_outlet, length_outlet, 0.0) +
                                     2.0 * Vecd(BW, BW, BW) +
                                     Vecd(extend_domain_length, extend_domain_length, 0.0));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real U_inlet = 1.0; //%
Real U_f = U_inlet; //*Characteristic velocity

Real U_max = 8.0 * U_inlet; //** An estimated value, this case one outlet and 8 inlets *

Real c_f = 10.0 * U_max;
Real rho0_f = 1.0; /**< Density. */

//Real Re = 3.2;
Real Re = 17.1;

Real Outlet_pressure = 0.0;

Real mu_f = rho0_f * U_f * D_hydraulic / Re;
Real Re_calculated = U_f * D_hydraulic * rho0_f / mu_f;

Real t_ref = 2.0; //% An estimated value

Real mixing_rate_interactive_radius = 1.4 * resolution_ref;
//----------------------------------------------------------------------
//	The open boundary setting.
//----------------------------------------------------------------------
/*
//% For G1
Vecd point_1(7.0, -8.0, 0.5);
Vecd point_2(10.6066, -0.707107, 0.5);
Vecd point_3(8.0, 7.0, 0.5);
Vecd point_4(0.707107, 10.6066, 0.5);
Vecd point_5(-7.0, 8.0, 0.5);
Vecd point_6(-10.6066, 0.707107, 0.5);
Vecd point_7(-8.0, -7.0, 0.5);
Vecd point_8(-0.707107, -10.6066, 0.5);
Vecd point_out(11.9373, 0.0, 3.5);
*/
/*
//% For G2
Vecd point_1(7.0, -5.3, 0.5);
Vecd point_2(8.69741, 1.20208, 0.5);
Vecd point_3(5.3, 7.0, 0.5);
Vecd point_4(-1.20208, 8.69741, 0.5);
Vecd point_5(-7.0, 5.3, 0.5);
Vecd point_6(-8.69741, -1.20208, 0.5);
Vecd point_7(-5.3, -7.0, 0.5);
Vecd point_8(1.20208, -8.69741, 0.5);
Vecd point_out(11.9373, 0.0, 3.5);
*/
/*
//% For G3
Vecd point_1(7.0, -5.5, 0.5);
Vecd point_2(8.83883, 1.06066, 0.5);
Vecd point_3(5.5, 7.0, 0.5);
Vecd point_4(-1.06066, 8.83883, 0.5);
Vecd point_5(-7.0, 5.5, 0.5);
Vecd point_6(-8.83883, -1.06066, 0.5);
Vecd point_7(-5.5, -7.0, 0.5);
Vecd point_8(1.06066, -8.83883, 0.5);
Vecd point_out(11.9373, 0.0, 3.5);
*/

/*
//% For G4
Vecd point_1(7.0, -6.0, 0.5);
Vecd point_2(9.19239, 0.707107, 0.5);
Vecd point_3(6.0, 7.0, 0.5);
Vecd point_4(-0.707107, 9.19239, 0.5);
Vecd point_5(-7.0, 6.0, 0.5);
Vecd point_6(-9.19239, -0.707107, 0.5);
Vecd point_7(-6.0, -7.0, 0.5);
Vecd point_8(0.707107, -9.19239, 0.5);
Vecd point_out(11.9373, 0.0, 3.5);
*/

/*
//% For G5
Vecd point_1(7.0, -6.0, 0.5);
Vecd point_2(9.19239, 0.707107, 0.5);
Vecd point_3(6.0, 7.0, 0.5);
Vecd point_4(-0.707107, 9.19239, 0.5);
Vecd point_5(-7.0, 6.0, 0.5);
Vecd point_6(-9.19239, -0.707107, 0.5);
Vecd point_7(-6.0, -7.0, 0.5);
Vecd point_8(0.707107, -9.19239, 0.5);
Vecd point_out(8.93725, 0.0, 3.5);
*/

/*
//% For G6
Vecd point_1(7.0, -6.0, 0.5);
Vecd point_2(9.19239, 0.707107, 0.5);
Vecd point_3(6.0, 7.0, 0.5);
Vecd point_4(-0.707107, 9.19239, 0.5);
Vecd point_5(-7.0, 6.0, 0.5);
Vecd point_6(-9.19239, -0.707107, 0.5);
Vecd point_7(-6.0, -7.0, 0.5);
Vecd point_8(0.707107, -9.19239, 0.5);
Vecd point_out(9.93725, 0.0, 3.5);
*/

/*
//% For G7
Vecd point_1(7.0, -6.0, 0.5);
Vecd point_2(9.19239, 0.707107, 0.5);
Vecd point_3(6.0, 7.0, 0.5);
Vecd point_4(-0.707107, 9.19239, 0.5);
Vecd point_5(-7.0, 6.0, 0.5);
Vecd point_6(-9.19239, -0.707107, 0.5);
Vecd point_7(-6.0, -7.0, 0.5);
Vecd point_8(0.707107, -9.19239, 0.5);
Vecd point_out(10.93725, 0.0, 3.5);
*/

/*
//% For G8
Vecd point_1(7.0, -6.0, 0.5);
Vecd point_2(9.19239, 0.707107, 0.5);
Vecd point_3(6.0, 7.0, 0.5);
Vecd point_4(-0.707107, 9.19239, 0.5);
Vecd point_5(-7.0, 6.0, 0.5);
Vecd point_6(-9.19239, -0.707107, 0.5);
Vecd point_7(-6.0, -7.0, 0.5);
Vecd point_8(0.707107, -9.19239, 0.5);
Vecd point_out(11.93725, 0.0, 3.5);
*/

/*
//% For G9
Vecd point_1(36.5, -31.0185, 2.0);
Vecd point_2(47.7428, 3.876, 2.0);
Vecd point_3(31.0185, 36.5, 2.0);
Vecd point_4(-3.876, 47.7428, 2.0);
Vecd point_5(-36.5, 31.0185, 2.0);
Vecd point_6(-47.7428, -3.876, 2.0);
Vecd point_7(-31.0185, -36.5, 2.0);
Vecd point_8(3.876, -47.7428, 2.0);
Vecd point_out(51.1977, 0.0, 98.0);
*/

/*
//% For G9 scale to m
Vecd point_1_temp(36.5, -31.0185, 2.0);
Vecd point_2_temp(47.7428, 3.876, 2.0);
Vecd point_3_temp(31.0185, 36.5, 2.0);
Vecd point_4_temp(-3.876, 47.7428, 2.0);
Vecd point_5_temp(-36.5, 31.0185, 2.0);
Vecd point_6_temp(-47.7428, -3.876, 2.0);
Vecd point_7_temp(-31.0185, -36.5, 2.0);
Vecd point_8_temp(3.876, -47.7428, 2.0);
Vecd point_out_temp(51.1977, 0.0, 98.0);

Real scale_temp = 1.0e-5;
Vecd point_1 = point_1_temp * scale_temp;
Vecd point_2 = point_2_temp * scale_temp;
Vecd point_3 = point_3_temp * scale_temp;
Vecd point_4 = point_4_temp * scale_temp;
Vecd point_5 = point_5_temp * scale_temp;
Vecd point_6 = point_6_temp * scale_temp;
Vecd point_7 = point_7_temp * scale_temp;
Vecd point_8 = point_8_temp * scale_temp;
Vecd point_out = point_out_temp * scale_temp;
*/

/*
//% For G11 scale to m
Vecd point_1_temp(36.5, -47.0185, 2.0);
Vecd point_2_temp(59.0565, -7.43771, 2.0);
Vecd point_3_temp(47.0185, 36.5, 2.0);
Vecd point_4_temp(7.43771, 59.0565, 2.0);
Vecd point_5_temp(-36.5, 47.0185, 2.0);
Vecd point_6_temp(-59.0565, 7.43771, 2.0);
Vecd point_7_temp(-47.0185, -36.5, 2.0);
Vecd point_8_temp(-7.43771, -59.0565, 2.0);
Vecd point_out_temp(51.1977, 0.0, 98.0);

Real scale_temp = 1.0e-5;
Vecd point_1 = point_1_temp * scale_temp;
Vecd point_2 = point_2_temp * scale_temp;
Vecd point_3 = point_3_temp * scale_temp;
Vecd point_4 = point_4_temp * scale_temp;
Vecd point_5 = point_5_temp * scale_temp;
Vecd point_6 = point_6_temp * scale_temp;
Vecd point_7 = point_7_temp * scale_temp;
Vecd point_8 = point_8_temp * scale_temp;
Vecd point_out = point_out_temp * scale_temp;
*/

///*
//% For G12
//Vecd point_1(36.5, -33.0185, 2.0);
//Vecd point_2(49.1570, 2.46178, 2.0);
//Vecd point_3(33.0185, 36.5, 2.0);
//Vecd point_4(-2.46178, 49.1570, 2.0);
//Vecd point_5(-36.5, 33.0185, 2.0);
//Vecd point_6(-49.1570, -2.46178, 2.0);
//Vecd point_7(-33.0185, -36.5, 2.0);
//Vecd point_8(2.46178, -49.1570, 2.0);
//Vecd point_out(51.1977, 0.0, 98.0);
//*/

///*
//% For G14
Vecd point_1(36.5, -42.0185, 2.0);
//Vecd point_2(55.5210, -3.90218, 2.0);
//Vecd point_3(42.0185, 36.5000, 2.0);
//Vecd point_4(3.90218, 55.5210, 2.0);
//Vecd point_5(-36.5000, 42.0185, 2.0);
//Vecd point_6(-55.5210, 3.90218, 2.0);
//Vecd point_7(-42.0185, -36.5000, 2.0);
//Vecd point_8(-3.90218, -55.5210, 2.0);
Vecd point_out(36.5, -27.0185, 2.0);


//% For G14 Initial color bounding point 
StdVec<Vecd> box_initial_bounding;
void getInitialBoundingBox()
{
    box_initial_bounding.push_back(point_1);
    //box_initial_bounding.push_back(point_2);
    //box_initial_bounding.push_back(point_3);
    //box_initial_bounding.push_back(point_4);
    //box_initial_bounding.push_back(point_5);
    //box_initial_bounding.push_back(point_6);
    //box_initial_bounding.push_back(point_7);
    //box_initial_bounding.push_back(point_8);
}
//*/


Vecd axis_vector_x(1.0, 0.0, 0.0);
Vecd axis_vector_z(0.0, 0.0, 1.0);
Vecd axis_vector_y(0.0, 1.0, 0.0);

//** L_inlet for X, axis-FLOW for Y, H_inlet for Z, */
Vecd inlet_buffer_halfsize = 0.5 * Vecd(L_inlet + 0.5 * resolution_ref, buffer_thickness, H_inlet);

Real inlet_1_rotation_angle = 0.0;
Rotation3d inlet_1_rotation(inlet_1_rotation_angle, axis_vector_z);
Eigen::AngleAxisd rotation_1(inlet_1_rotation_angle, axis_vector_z);
Vecd inlet_1_flow_unit_vector = rotation_1 * axis_vector_y;
Vecd inlet_1_buffer_translation = point_1 + 0.5 * buffer_thickness * inlet_1_flow_unit_vector;
Vecd inlet_1_sub_buffer_translation = inlet_1_buffer_translation - buffer_thickness * inlet_1_flow_unit_vector;
//
//Real inlet_2_rotation_angle = 45.0 * M_PI / 180.0;
//Rotation3d inlet_2_rotation(inlet_2_rotation_angle, axis_vector_z);
//Eigen::AngleAxisd rotation_2(inlet_2_rotation_angle, axis_vector_z);
//Vecd inlet_2_flow_unit_vector = rotation_2 * axis_vector_y;
//Vecd inlet_2_buffer_translation = point_2 + 0.5 * buffer_thickness * inlet_2_flow_unit_vector;
//Vecd inlet_2_sub_buffer_translation = inlet_2_buffer_translation - buffer_thickness * inlet_2_flow_unit_vector;
//
//Real inlet_3_rotation_angle = 45.0 * 2.0 * M_PI / 180.0;
//Rotation3d inlet_3_rotation(inlet_3_rotation_angle, axis_vector_z);
//Eigen::AngleAxisd rotation_3(inlet_3_rotation_angle, axis_vector_z);
//Vecd inlet_3_flow_unit_vector = rotation_3 * axis_vector_y;
//Vecd inlet_3_buffer_translation = point_3 + 0.5 * buffer_thickness * inlet_3_flow_unit_vector;
//Vecd inlet_3_sub_buffer_translation = inlet_3_buffer_translation - buffer_thickness * inlet_3_flow_unit_vector;
//
//Real inlet_4_rotation_angle = 45.0 * 3.0 * M_PI / 180.0;
//Rotation3d inlet_4_rotation(inlet_4_rotation_angle, axis_vector_z);
//Eigen::AngleAxisd rotation_4(inlet_4_rotation_angle, axis_vector_z);
//Vecd inlet_4_flow_unit_vector = rotation_4 * axis_vector_y;
//Vecd inlet_4_buffer_translation = point_4 + 0.5 * buffer_thickness * inlet_4_flow_unit_vector;
//Vecd inlet_4_sub_buffer_translation = inlet_4_buffer_translation - buffer_thickness * inlet_4_flow_unit_vector;
//
//Real inlet_5_rotation_angle = 45.0 * 4.0 * M_PI / 180.0;
//Rotation3d inlet_5_rotation(inlet_5_rotation_angle, axis_vector_z);
//Eigen::AngleAxisd rotation_5(inlet_5_rotation_angle, axis_vector_z);
//Vecd inlet_5_flow_unit_vector = rotation_5 * axis_vector_y;
//Vecd inlet_5_buffer_translation = point_5 + 0.5 * buffer_thickness * inlet_5_flow_unit_vector;
//Vecd inlet_5_sub_buffer_translation = inlet_5_buffer_translation - buffer_thickness * inlet_5_flow_unit_vector;
//
//Real inlet_6_rotation_angle = 45.0 * 5.0 * M_PI / 180.0;
//Rotation3d inlet_6_rotation(inlet_6_rotation_angle, axis_vector_z);
//Eigen::AngleAxisd rotation_6(inlet_6_rotation_angle, axis_vector_z);
//Vecd inlet_6_flow_unit_vector = rotation_6 * axis_vector_y;
//Vecd inlet_6_buffer_translation = point_6 + 0.5 * buffer_thickness * inlet_6_flow_unit_vector;
//Vecd inlet_6_sub_buffer_translation = inlet_6_buffer_translation - buffer_thickness * inlet_6_flow_unit_vector;
//
//Real inlet_7_rotation_angle = 45.0 * 6.0 * M_PI / 180.0;
//Rotation3d inlet_7_rotation(inlet_7_rotation_angle, axis_vector_z);
//Eigen::AngleAxisd rotation_7(inlet_7_rotation_angle, axis_vector_z);
//Vecd inlet_7_flow_unit_vector = rotation_7 * axis_vector_y;
//Vecd inlet_7_buffer_translation = point_7 + 0.5 * buffer_thickness * inlet_7_flow_unit_vector;
//Vecd inlet_7_sub_buffer_translation = inlet_7_buffer_translation - buffer_thickness * inlet_7_flow_unit_vector;
//
//Real inlet_8_rotation_angle = 45.0 * 7.0 * M_PI / 180.0;
//Rotation3d inlet_8_rotation(inlet_8_rotation_angle, axis_vector_z);
//Eigen::AngleAxisd rotation_8(inlet_8_rotation_angle, axis_vector_z);
//Vecd inlet_8_flow_unit_vector = rotation_8 * axis_vector_y;
//Vecd inlet_8_buffer_translation = point_8 + 0.5 * buffer_thickness * inlet_8_flow_unit_vector;
//Vecd inlet_8_sub_buffer_translation = inlet_8_buffer_translation - buffer_thickness * inlet_8_flow_unit_vector;

//** axis-FLOW for X, L_outlet for Y, H_outlet for Z,  */
Vecd outlet_buffer_halfsize = 0.5 * Vecd(L_inlet + 0.5 * resolution_ref, buffer_thickness, H_inlet);

Real outlet_rotation_angle = 0.0 * M_PI / 180.0;
Real outlet_rotation_angle_reverse = M_PI + outlet_rotation_angle; //%% In this code, the outlet is defined by the inflow direction

Rotation3d outlet_rotation(outlet_rotation_angle, axis_vector_z);
Rotation3d outlet_rotation_reverse(outlet_rotation_angle_reverse, axis_vector_z);

Eigen::AngleAxisd rotation_outlet(outlet_rotation_angle, axis_vector_z);
Vecd outlet_flow_unit_vector = rotation_outlet * axis_vector_y;
Vecd outlet_buffer_translation = point_out - 0.5 * buffer_thickness * outlet_flow_unit_vector;
Vecd outlet_sub_buffer_translation = outlet_buffer_translation + buffer_thickness * outlet_flow_unit_vector;

// Vecd left_buffer_translation = Vecd(0.0, 0.0, 0.0); //** STL relevant */
// Vecd right_buffer_halfsize = 0.5 * Vecd(L_outlet, H_outlet, buffer_thickness);
// Vecd right_buffer_translation = Vecd(0.0, 0.0, 0.0); //** STL relevant */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::string stl_fluid_path = "./input/g14-return-straight-0p19.stl";
Real scale_factor_fluid = 1.0;
Vecd translation_stl_fluid(0.0, 0.0, 0.0);
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stl_fluid_path, translation_stl_fluid, scale_factor_fluid);
    }
};

/** Set the file path to the stl file. */
std::string stl_structure_path = "./input/g14-return-straight-0p19.stl"; //% This also denote which file we use
Real scale_factor = 1.0;
Vecd translation_stl(0.0, 0.0, 0.0);
class WallBoundaryFromSTL : public ComplexShape
{
  public:
    explicit WallBoundaryFromSTL(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<ExtrudeShape<TriangleMeshShapeSTL>>(BW, stl_structure_path, translation_stl, scale_factor);
        subtract<TriangleMeshShapeSTL>(stl_structure_path, translation_stl, scale_factor);
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_;
    AlignedBox &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_inlet),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = current_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref)) : u_ref_;
        target_velocity = u_ave * axis_vector_y;
        return target_velocity;
    }
};
struct RightOutflowPressure
{
    template <class BoundaryConditionType>
    RightOutflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
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

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};
//----------------------------------------------------------------------
// Observation with offset model.
//----------------------------------------------------------------------
//** For getting centerline velocity *
// namespace observe_centerline
// {
// constexpr const char *namespace_prefix = "centerline";
// Vecd pos_observe_start = point_O;
// Real sparsity_ratio = 5.0;
// Real length_observing_line = DL;
// Vecd unit_direction_observe = Vecd(0.0, 0.0, 1.0);
// Real observer_offset_distance = 0.0;

// int num_observer_points = std::round(length_observing_line / resolution_ref / sparsity_ratio); //**Every particle is regarded as a cell monitor*
// Real observe_spacing = length_observing_line / num_observer_points;

// StdVec<Vecd> observation_location;
// void get_observation_locations()
// {
//     for (int i = 0; i < num_observer_points; ++i)
//     {
//         Vecd pos_observer_i = pos_observe_start + i * observe_spacing * unit_direction_observe;
//         if (i == 0)
//         {
//             pos_observer_i -= observer_offset_distance * unit_direction_observe;
//         }
//         if (i == num_observer_points - 1)
//         {
//             pos_observer_i += observer_offset_distance * unit_direction_observe;
//         }
//         observation_location.push_back(pos_observer_i);
//     }
// }
// void output_observer_theoretical_pos_on_line()
// {
//     std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_theoretical_pos_on_line.dat";
//     std::ofstream outfile(filename);
//     if (!outfile.is_open())
//     {
//         std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
//         return;
//     }
//     for (int i = 0; i < num_observer_points; ++i)
//     {
//         outfile << observation_location[i].dot(unit_direction_observe) << "\n";
//     }
//     outfile.close();
// }
// } // namespace observe_centerline

// //** For getting cross-section velocity *
// namespace observe_cross_sections
// {
// constexpr const char *namespace_prefix = "cross_sections";
// const int number_observe_line = 5;
// Real observer_offset_distance = 2.0 * resolution_ref;
// Vecd unit_direction_observe(0.0, 1.0, 0.0);
// // ** Determine the observing start point. *
// Real observe_start_z[number_observe_line] = {
//     0.0 * DL + observer_offset_distance,
//     0.25 * DL,
//     0.50 * DL,
//     0.75 * DL,
//     0.99 * DL - observer_offset_distance};

// Real observe_start_x[number_observe_line] = {
//     0.0,
//     0.0,
//     0.0,
//     0.0,
//     0.0};

// Real observe_start_y[number_observe_line] = {
//     0.5 * resolution_ref - Radius_inlet,
//     0.5 * resolution_ref - Radius_inlet,
//     0.5 * resolution_ref - Radius_inlet,
//     0.5 * resolution_ref - Radius_inlet,
//     0.5 * resolution_ref - Radius_inlet};

// // ** Determine the length of the observing line and other information. *
// Real observe_line_length[number_observe_line] = {0.0};
// int num_observer_points[number_observe_line] = {0};

// void getObservingLineLengthAndEndPoints()
// {
//     for (int i = 0; i < number_observe_line; ++i)
//     {
//         observe_line_length[i] = DH;
//         num_observer_points[i] = std::round(observe_line_length[i] / resolution_ref);
//     }
// }

// StdVec<Vecd> observation_locations;
// StdVec<Vecd> observation_theoretical_locations;
// void getPositionsOfMultipleObserveLines()
// {
//     getObservingLineLengthAndEndPoints();
//     for (int k = 0; k < number_observe_line; ++k)
//     {
//         Vecd pos_observe_start(observe_start_x[k], observe_start_y[k], observe_start_z[k]);
//         int num_observer_point = num_observer_points[k];
//         Real observe_spacing = observe_line_length[k] / num_observer_point;
//         for (int i = 0; i < num_observer_point; ++i)
//         {
//             Real offset = 0.0;
//             offset = (i == 0 ? -observer_offset_distance : (i == num_observer_point - 1 ? observer_offset_distance : 0.0));
//             Vecd pos_observer_i = pos_observe_start + (i * observe_spacing + offset) * unit_direction_observe;
//             Vecd pos_observer_i_no_offset = pos_observe_start + i * observe_spacing * unit_direction_observe;
//             observation_locations.push_back(pos_observer_i);
//             observation_theoretical_locations.push_back(pos_observer_i_no_offset);
//         }
//     }
// }
// void output_observe_positions()
// {
//     std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_observer_positions.dat";
//     std::ofstream outfile(filename);
//     if (!outfile.is_open())
//     {
//         std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
//         return;
//     }
//     for (const Vecd &position : observation_locations)
//     {
//         for (int i = 0; i < position.size(); ++i)
//         {
//             outfile << position[i] << " ";
//         }
//         outfile << "\n";
//     }
//     outfile.close();
// }
// void output_observer_theoretical_pos_on_line()
// {
//     std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_theoretical_pos_on_line.dat";
//     std::ofstream outfile(filename);
//     if (!outfile.is_open())
//     {
//         std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
//         return;
//     }
//     for (int j = 0; j < number_observe_line; ++j)
//     {
//         for (int i = 0; i < num_observer_points[j]; ++i)
//         {
//             outfile << observation_theoretical_locations[i].dot(unit_direction_observe) << "\n";
//         }
//     }
//     outfile.close();
// }
// void output_number_observe_points_on_lines()
// {
//     std::string filename = "../bin/output/" + std::string(namespace_prefix) + "_observer_num_points_on_lines.dat";
//     std::ofstream outfile(filename);
//     if (!outfile.is_open())
//     {
//         std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
//         return;
//     }
//     for (const int &number : num_observer_points)
//     {
//         outfile << number << "\n";
//     }
//     outfile.close();
// }
// } // namespace observe_cross_sections

//** For regression test *
StdVec<Vecd> observer_location_center_point = {point_O + Vecd(0.0, 0.0, 0.5 * H_total)};