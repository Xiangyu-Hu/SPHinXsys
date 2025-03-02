/**
 * @file sloshing_hanging_baffle.h
 * @The benchmark example of a rolling tank with the 
   elastic baffle hanging on the top.
 * @author Bo Zhang and Xiangyu Hu
 */

#ifndef SLOSHING_HANGING_BAFFLE_H
#define SLOSHING_HANGING_BAFFLE_H

#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.608;                         /**< Channel length. */
Real DH = 0.3445;                        /**< Channel height. */
Real Water_height = 0.1148;              /**< Water block height. */
Real Baffle_width = 0.004;               /**< Hanging baggle width. */
Real particle_spacing_ref = 0.002;       /**< Global reference resolution. */
Real BW = 4.0 * particle_spacing_ref;    /**< Extending width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW - 0.5 * DL, -BW), Vec2d(0.5 * DL + BW, DH + BW));
/**< Original point is in left bottom. */
//----------------------------------------------------------------------
//	Define the corner point of water block geometry.
//----------------------------------------------------------------------
Vec2d DamP_lb(-DL / 2, 0.0);               /**< Left bottom. */
Vec2d DamP_lt(-DL / 2, Water_height);      /**< Left top. */
Vec2d DamP_rt(DL / 2, Water_height);       /**< Right top. */
Vec2d DamP_rb(DL / 2, 0.0);                /**< Right bottom. */
//----------------------------------------------------------------------
//	Define the corner point of gate geometry.
//----------------------------------------------------------------------
Vec2d Baffle_lb(-0.5 * Baffle_width, -BW);
Vec2d Baffle_lt(-0.5 * Baffle_width, Water_height);
Vec2d Baffle_rt(0.5 * Baffle_width, Water_height);
Vec2d Baffle_rb(0.5 * Baffle_width, -BW);
//----------------------------------------------------------------------
//	Define the geometry for gate constrain.
//----------------------------------------------------------------------
Vec2d Constrain_lb(-0.5 * Baffle_width, -BW);
Vec2d Constrain_lt(-0.5 * Baffle_width, 0);
Vec2d Constrain_rt(0.5 * Baffle_width, 0);
Vec2d Constrain_rb(0.5 * Baffle_width, -BW);
// observer location
StdVec<Vecd> observation_location = { Vecd(0.0, Water_height) };
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 917.0;   /**< Reference density of fluid. */
Real gravity_g = 9.81; /**< Value of gravity. */
Real U_ref = 2.0 * sqrt(Water_height * gravity_g);
;                                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;              /**< Reference sound speed. */
Real mu_f = 0.00005 * rho0_f;         /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Material properties of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 1100.0; /**< Reference solid density. */
Real poisson = 0.49;  /**< Poisson ratio. */
Real Ae = 6e6;        /**< Normalized Youngs Modulus : 4.0 MPa. */
Real Youngs_modulus = Ae;
//----------------------------------------------------------------------
//	Define the baffle body shape
//----------------------------------------------------------------------
/** create a baffle shape */
std::vector<Vecd> createBaffleShape()
{
    std::vector<Vecd> baffle_shape;
    baffle_shape.push_back(Baffle_lb);
    baffle_shape.push_back(Baffle_lt);
    baffle_shape.push_back(Baffle_rt);
    baffle_shape.push_back(Baffle_rb);
    baffle_shape.push_back(Baffle_lb);

    return baffle_shape;
}
class Baffle : public MultiPolygonShape
{
public:
    explicit Baffle(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createBaffleShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(DamP_lb);
    water_block_shape.push_back(DamP_lt);
    water_block_shape.push_back(DamP_rt);
    water_block_shape.push_back(DamP_rb);
    water_block_shape.push_back(DamP_lb);

    return water_block_shape;
}
class WaterBlock : public MultiPolygonShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createBaffleShape(), ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	create baffle constrain shape
//----------------------------------------------------------------------
std::vector<Vecd> createBaffleConstrainShape()
{
    // geometry
    std::vector<Vecd> gate_constraint_shape;
    gate_constraint_shape.push_back(Constrain_lb);
    gate_constraint_shape.push_back(Constrain_lt);
    gate_constraint_shape.push_back(Constrain_rt);
    gate_constraint_shape.push_back(Constrain_rb);
    gate_constraint_shape.push_back(Constrain_lb);

    return gate_constraint_shape;
};

MultiPolygon BaffleConstraint()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(createBaffleConstrainShape(), ShapeBooleanOps::add);
    return multi_polygon;
};
//----------------------------------------------------------------------
//	create wall boundary shape
//----------------------------------------------------------------------
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-BW - 0.5 * DL, -BW));
    outer_wall_shape.push_back(Vecd(-BW - 0.5 * DL, DH + BW));
    outer_wall_shape.push_back(Vecd(BW + 0.5 * DL, DH + BW));
    outer_wall_shape.push_back(Vecd(BW + 0.5 * DL, -BW));
    outer_wall_shape.push_back(Vecd(-BW - 0.5 * DL, -BW));

    return outer_wall_shape;
}
/** create inner wall shape */
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(-0.5 * DL, 0.0));
    inner_wall_shape.push_back(Vecd(-0.5 * DL, DH));
    inner_wall_shape.push_back(Vecd(0.5 * DL, DH));
    inner_wall_shape.push_back(Vecd(0.5 * DL, 0.0));
    inner_wall_shape.push_back(Vecd(-0.5 * DL, 0.0));

    return inner_wall_shape;
}

class WallBoundary : public MultiPolygonShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createBaffleConstrainShape(), ShapeBooleanOps::sub);
    }
}; 

Real omega = 2.0 * PI * 0.82576; // period of time;
Real Theta0 = 4.0 * PI / 180;   // maximum rotating angle;

class VariableGravity : public Gravity
{
public:
    VariableGravity(Vecd gravity_vector) : Gravity(gravity_vector) {};
    Vecd InducedAcceleration(const Vecd& position, Real physical_time) const
    {
        Vecd reference_acceleration_ = Vecd::Zero();

        Real Theta = Theta0 * sin(omega * (physical_time - 0.3));
        Real ThetaV = Theta0 * omega * cos(omega * (physical_time - 0.3));
        
        Real alpha = std::atan2(position[1], position[0]);
        Real distance = std::sqrt(pow(position[0], 2) + pow(position[1], 2));

        if (physical_time <= 0.3) 
        {
            reference_acceleration_[0] = 0.0;
            reference_acceleration_[1] = -gravity_g;
        }
        else {
            // 重力分量
            reference_acceleration_[0] = -gravity_g * sin(Theta);
            reference_acceleration_[1] = -gravity_g * cos(Theta);

            // 离心力项
            reference_acceleration_[0] -= ThetaV * ThetaV * position[0];
            reference_acceleration_[1] += ThetaV * ThetaV * position[1];

            // 科氏力项
            Real Vx = -ThetaV * position[1];
            Real Vy = ThetaV * position[0];

            reference_acceleration_[0] += 2 * ThetaV * Vy; // 2ω_z v'_y
            reference_acceleration_[1] -= 2 * ThetaV * Vx; // -2ω_z v'_x
        }

        return reference_acceleration_;
    }
};
// namespace SPH
#endif // SLOSHING_HANGING_BAFFLE_H
