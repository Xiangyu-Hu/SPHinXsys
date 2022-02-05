/**
* @file 	2d_eulerian_flow_around_cylinder.h
* @brief 	This is the case file for the test of flow passing by a cylinder
			in eulerian framework.
* @details  We consider a flow passing by a cylinder in 2D.
* @author 	Zhentong Wang,Xiangyu Hu and Chi Zhang
*/

#ifndef _2D_EULERIAN_FLOW_AROUND_CYLINDER_H
#define _2D_EULERIAN_FLOW_AROUND_CYLINDER_H

#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters.
 */
Real DL = 15.0; 					            /**< Channel length. */
Real DH = 10.0; 						        /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; 	            /**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 2.0;			/**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0;			/**< Sponge region to impose inflow condition. */
Vec2d insert_circle_center(4.0, 5.0);		    /**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;			    /**< Radius of the cylinder. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));

/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		/**< Density. */
Real U_f = 1.0;			/**< freestream velocity. */
Real c_f = 10.0 * U_f;	/**< Speed of sound. */
Real Re = 100.0;		/**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re;	/**< Dynamics viscosity. */

/**
* @brief define geometry of SPH bodies
*/
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

	return water_block_shape;
}

/**
 * @brief Define case dependent bodies material, constraint and boundary conditions.
 */
 /** Fluid body definition */
class WaterBlock : public EulerianFluidBody
{
public:
	WaterBlock(SPHSystem& system, std::string body_name)
		: EulerianFluidBody(system, body_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};

/** Definition of the cylinder. */
class Cylinder : public SolidBody
{
public:
	Cylinder(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 2.0))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};

/** far field boundary condition and set up parameters. */
class FarFieldBoundary : public eulerian_weakly_compressible_fluid_dynamics::NonReflectiveBoundaryVariableCorrection
{
public:
	FarFieldBoundary(BaseBodyRelationInner &inner_relation) :
		eulerian_weakly_compressible_fluid_dynamics::NonReflectiveBoundaryVariableCorrection(inner_relation)
	{
		rho_farfield_ = rho0_f;
		vel_farfield_ = Vecd(U_f, 0.0);
		sound_speed_ = c_f;
	};
	virtual ~FarFieldBoundary() {};
};
#endif //_2D_FLOW_AROUND_CYLINDER_H
