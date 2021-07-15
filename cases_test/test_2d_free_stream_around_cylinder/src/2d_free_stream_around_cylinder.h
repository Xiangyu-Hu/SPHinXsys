/**
* @file 	freestream_flow_around_cylinder_case.h
* @brief 	This is the case file for the test of free-stream flow.
* @details  We consider a flow pass the cylinder with freestream boundary condition in 2D.
* @author 	Xiangyu Hu, Shuoguo Zhang
*/

#pragma once

#include "sphinxsys.h"

using namespace SPH;
/**
 * @brief Basic geometry parameters.
 */
Real DL = 30.0; 					/**< Channel length. */
Real DH = 16.0; 				    /**< Channel height. */
Real particle_spacing_ref = 0.2; 	/**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0;	/**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;           /**< Sponge region to impose injection. */
Vec2d insert_circle_center(10.0, 0.5 * DH);		/**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;			    /**< Radius of the cylinder. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
/** Prescribed fluid body domain bounds*/
BoundingBox fluid_body_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		/**< Density. */
Real U_f = 1.0;			/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;	/**< Speed of sound. */
Real Re = 100.0;		/**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re;	/**< Dynamics viscosity. */
/**
* @brief define fluid
*/
/** create a water block shape */
std::vector<Vecd> CreatWaterBlockShape()
{
	/** Geomtry definition. */
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, DH));
	water_block_shape.push_back(Vecd(DL, DH));
	water_block_shape.push_back(Vecd(DL, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

	return water_block_shape;
}
/** Fluid body definition */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& system, string body_name)
		: FluidBody(system, body_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> water_block_shape = CreatWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);
		body_shape_->addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
	}
};
/** fluid material properties. */
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() :WeaklyCompressibleFluid()
	{
		rho0_ = rho0_f;
		c0_ = c_f;
		mu_ = mu_f;

		/** supplementary material paramters derived from basic parameters. */
		assignDerivedMaterialParameters();
	}
};
/**
* @brief define solid
*/
/**  create a cylinder shape and define solid body. */
class Cylinder : public SolidBody
{
public:
	Cylinder(SPHSystem& system, string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(1.15, 1))
	{
		/** Geomtry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
/**
* @brief define the emitter buffer
*/
/** create the emitter buffer shape . */
std::vector<Vecd> CreatEmitterBufferShape()
{
	/** Geomtry definition. */
	std::vector<Vecd> emitter_buffer_shape;
	emitter_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
	emitter_buffer_shape.push_back(Vecd(-DL_sponge, DH));
	emitter_buffer_shape.push_back(Vecd(-0.0, DH));
	emitter_buffer_shape.push_back(Vecd(-0.0, 0.0));
	emitter_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

	return emitter_buffer_shape;
}
/** define emitter buffer body*/
class EmitterBuffer : public BodyPartByParticle
{
public:
	EmitterBuffer(FluidBody* fluid_body, string constrained_region_name)
		: BodyPartByParticle(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> emitter_buffer_shape = CreatEmitterBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(emitter_buffer_shape, ShapeBooleanOps::add);

		//tag the constrained particle
		tagBodyPart();
	}
};
/** define emitter inflow boundary conditon*/
class EmitterInflowCondition : public fluid_dynamics::InletOutletInflowCondition
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterInflowCondition(FluidBody* fluid_body,
		BodyPartByParticle* constrained_region)
		: InletOutletInflowCondition(fluid_body, constrained_region)
	{
		u_ave_ = 0.0;
		u_ref_ = U_f;
		t_ref_ = 2.0;
		SetInflowParameters();
	}

	Vecd getTargetVelocity(Vecd& position, Vecd& velocity)
	{
		Real u = velocity[0];
		Real v = velocity[1];

		if (position[0] < 0.0)
		{
			u = u_ave_;
			v = 0.0;
		}
		return Vecd(u, v);
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
	}

	void SetInflowParameters() override
	{
		inflow_pressure_ = 0.0;
	}
};
/**
* @brief define the injection buffer
*/
/** create the injection buffer shape. */
std::vector<Vecd> CreatInjectionBufferShape()
{
	/** Geomtry definition. */
	std::vector<Vecd> injection_buffer_shape;
	injection_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
	injection_buffer_shape.push_back(Vecd(-DL_sponge, DH));
	injection_buffer_shape.push_back(Vecd(-DL_sponge + BW, DH));
	injection_buffer_shape.push_back(Vecd(-DL_sponge + BW, 0.0));
	injection_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

	return injection_buffer_shape;
}
/** define injection buffer body*/
class InjectionBuffer : public BodyPartByParticle
{
public:
	InjectionBuffer(FluidBody* fluid_body, string constrained_region_name)
		: BodyPartByParticle(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> injection_buffer_shape = CreatInjectionBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(injection_buffer_shape, ShapeBooleanOps::add);

		//tag the constrained particle
		tagBodyPart();
	}
};
/**
 * define time dependent acceleration in x-direction
 */
class TimeDependentAcceleration : public Gravity
{
	Real du_ave_dt_, u_ref_, t_ref_;
public:
	TimeDependentAcceleration(Vecd gravity_vector)
		: Gravity(gravity_vector)
	{
		t_ref_ = 2.0;
		u_ref_ = U_f;
		du_ave_dt_ = 0.0;
	}
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		Real run_time_ = GlobalStaticVariables::physical_time_;
		du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);

		return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
	}
};
/**
 * define fluid observer
 */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem& system, string body_name)
		: FictitiousBody(system, body_name)
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