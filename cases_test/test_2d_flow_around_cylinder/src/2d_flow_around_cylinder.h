/**
* @file 	2d_cylinder_flow.h
* @brief 	This is the case file for the test of flow passing by a cylinder.
* @details  We consider a flow passing by a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
* @version 0.1
*/

#pragma once

#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters.
 */
Real DL = 15.0; 					            /**< Channel length. */
Real DH = 10.0; 						        /**< Channel height. */
Real particle_spacing_ref = 0.2; 	            /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 10.0;	/**< Sponge region to impose inflow condition. */
Real DH_sponge = particle_spacing_ref * 2.0;					/**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0; 		    /**< Boundary width, determined by specific layer of boundary particles. */
Vec2d insert_circle_center(4.0, 5.0);		    /**< Location of the cylinder center. */
Real insert_circle_radius = 0.75;			    /**< Radius of the cylinder. */

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
std::vector<Point> CreatWaterBlockShape()
{
	//geometry
	std::vector<Point> water_block_shape;
	water_block_shape.push_back(Point(-DL_sponge, -DH_sponge));
	water_block_shape.push_back(Point(-DL_sponge, DH + DH_sponge));
	water_block_shape.push_back(Point(DL, DH + DH_sponge));
	water_block_shape.push_back(Point(DL, -DH_sponge));
	water_block_shape.push_back(Point(-DL_sponge, -DH_sponge));

	return water_block_shape;
}

/** create a water block buffer shape. */
std::vector<Point> CreatBufferShape()
{
	std::vector<Point> buffer_shape;
	buffer_shape.push_back(Point(-DL_sponge, -DH_sponge));
	buffer_shape.push_back(Point(-DL_sponge, DH + DH_sponge));
	buffer_shape.push_back(Point(DL, DH + DH_sponge));
	buffer_shape.push_back(Point(DL, DH));
	buffer_shape.push_back(Point(0.0, DH));
	buffer_shape.push_back(Point(0.0, 0.0));
	buffer_shape.push_back(Point(DL, 0.0));
	buffer_shape.push_back(Point(DL, -DH_sponge));
	buffer_shape.push_back(Point(-DL_sponge, -DH_sponge));

	return buffer_shape;
}
/**
 * @brief Define case dependent bodies material, constraint and boundary conditions.
 */
 /** Fluid body definition */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& system, string body_name, int refinement_level)
		: FluidBody(system, body_name, refinement_level)
	{
		/** Geomtry definition. */
		std::vector<Point> water_block_shape = CreatWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);
		body_shape_->addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
	}
};
/** Case-dependent material properties. */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		rho_0_ = rho0_f;
		c_0_ = c_f;
		mu_ = mu_f;
		/** supplementary material paramters derived from basic parameters. */
		assignDerivedMaterialParameters();
	}
};
/** Definition of the cylinder. */
class Cylinder : public SolidBody
{
public:
	Cylinder(SPHSystem& system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geomtry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
/** inflow buffer */
class FreeStreamBuffer : public BodyPartByCell
{
public:
	FreeStreamBuffer(FluidBody* fluid_body, string constrained_region_name)
		: BodyPartByCell(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Point> inflow_buffer_shape = CreatBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);

		//tag the constrained particle
		tagBodyPart();
	}
};
/** Case dependent inflow boundary condition. */
class FreeStreamCondition : public fluid_dynamics::FlowRelaxationBuffer
{
	Real u_ave_, u_ref_, t_ref;
public:
	FreeStreamCondition(FluidBody* fluid_body,
		BodyPartByCell* constrained_region)
		: fluid_dynamics::FlowRelaxationBuffer(fluid_body, constrained_region)
	{
		u_ave_ = 0.0;
		u_ref_ = U_f;
		t_ref = 2.0;
	}
	Vecd getTargetVelocity(Vecd& position, Vecd& velocity)
	{
		return Vecd(u_ave_, 0.0);
	}
	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;
	}
};
/** an observer body to measure the flow profile */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem& system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
	{
		/** the measureing particles */
		Vec2d point_coordinate_1(3.0, 5.0);
		Vec2d point_coordinate_2(4.0, 5.0);
		Vec2d point_coordinate_3(5.0, 5.0);

		body_input_points_volumes_.push_back(make_pair(point_coordinate_1, 0.0));
		body_input_points_volumes_.push_back(make_pair(point_coordinate_2, 0.0));
		body_input_points_volumes_.push_back(make_pair(point_coordinate_3, 0.0));
	}
};
