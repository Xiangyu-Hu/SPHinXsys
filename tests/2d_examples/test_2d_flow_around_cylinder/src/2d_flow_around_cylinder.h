/**
* @file 	2d_cylinder_flow.h
* @brief 	This is the case file for the test of flow passing by a cylinder.
* @details  We consider a flow passing by a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
*/

#ifndef _2D_FLOW_AROUND_CYLINDER_H
#define _2D_FLOW_AROUND_CYLINDER_H

#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters.
 */
Real DL = 15.0; 					            /**< Channel length. */
Real DH = 10.0; 						        /**< Channel height. */
Real resolution_ref = 0.2; 	            		/**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 10.0;			/**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0;			/**< Sponge region to impose inflow condition. */
Vec2d insert_circle_center(4.0, 5.0);		    /**< Location of the cylinder center. */
Real insert_circle_radius = 0.75;			    /**< Radius of the cylinder. */
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

/** create a water block buffer shape. */
std::vector<Vecd> CreatBufferShape()
{
	std::vector<Vecd> buffer_shape;
	buffer_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
	buffer_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
	buffer_shape.push_back(Vecd(DL, DH + DH_sponge));
	buffer_shape.push_back(Vecd(DL, DH));
	buffer_shape.push_back(Vecd(0.0, DH));
	buffer_shape.push_back(Vecd(0.0, 0.0));
	buffer_shape.push_back(Vecd(DL, 0.0));
	buffer_shape.push_back(Vecd(DL, -DH_sponge));
	buffer_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

	return buffer_shape;
}
/**
 * @brief Define case dependent bodies material, constraint and boundary conditions.
 */
 /** Fluid body definition */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& system, std::string body_name)
		: FluidBody(system, body_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> water_block_shape = createWaterBlockShape();
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
		rho0_ = rho0_f;
		c0_ = c_f;
		mu_ = mu_f;
		/** supplementary material paramters derived from basic parameters. */
		assignDerivedMaterialParameters();
	}
};

class ParameterizedWaterMaterial : public BaseParameterization<WaterMaterial>
{
public:
	ParameterizedWaterMaterial(ParameterizationIO& parameterization_io) : 
		BaseParameterization<WaterMaterial>(parameterization_io)
	{
		getAParameter("WaterMaterial", "Viscosity", mu_);
	}
};

/** Definition of the cylinder. */
class Cylinder : public SolidBody
{
public:
	Cylinder(SPHSystem& system, std::string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(1.15, 2.0))
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
	FreeStreamBuffer(FluidBody* fluid_body, std::string constrained_region_name)
		: BodyPartByCell(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> inflow_buffer_shape = CreatBufferShape();
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
	FluidObserver(SPHSystem& system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		/** the measureing particles */
		Vec2d point_coordinate_1(3.0, 5.0);
		Vec2d point_coordinate_2(4.0, 5.0);
		Vec2d point_coordinate_3(5.0, 5.0);

		body_input_points_volumes_.push_back(std::make_pair(point_coordinate_1, 0.0));
		body_input_points_volumes_.push_back(std::make_pair(point_coordinate_2, 0.0));
		body_input_points_volumes_.push_back(std::make_pair(point_coordinate_3, 0.0));
	}
};
#endif //_2D_FLOW_AROUND_CYLINDER_H