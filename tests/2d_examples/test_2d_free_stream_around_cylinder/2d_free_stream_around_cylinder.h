/**
* @file 	freestream_flow_around_cylinder_case.h
* @brief 	This is the case file for the test of free-stream boundary condition.
* @details  We consider a flow pass the cylinder with freestream boundary condition in 2D.
* @author 	Xiangyu Hu, Shuoguo Zhang
*/

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 30.0;								  /**< Channel length. */
Real DH = 16.0;								  /**< Channel height. */
Real particle_spacing_ref = 0.2;			  /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;		  /**< Sponge region to impose injection. */
Vec2d insert_circle_center(10.0, 0.5 * DH);	  /**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;			  /**< Radius of the cylinder. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
/** Prescribed fluid body domain bounds*/
BoundingBox fluid_body_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;											/**< Density. */
Real U_f = 1.0;												/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;										/**< Speed of sound. */
Real Re = 100.0;											/**< Reynolds number. */
Real mu = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
//	water block shape
std::vector<Vecd> water_block_shape {
Vecd(-DL_sponge, 0.0),Vecd(-DL_sponge, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(-DL_sponge, 0.0)};
//----------------------------------------------------------------------
//	Define case dependent SPH bodies.
//----------------------------------------------------------------------
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/**  define solid body. */
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
//----------------------------------------------------------------------
//	Define case dependent SPH body parts.
//----------------------------------------------------------------------
/** create the emitter shape. */
MultiPolygon creatEmitterShape()
{
	std::vector<Vecd> emmiter_shape{
		Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(-DL_sponge + BW, DH), Vecd(-DL_sponge + BW, 0.0), Vecd(-DL_sponge, 0.0)};

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(emmiter_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
/** create the emitter buffer shape . */
MultiPolygon createEmitterBufferShape()
{
	std::vector<Vecd> emitter_buffer_shape{
		Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(0.0, DH), Vecd(0.0, 0.0), Vecd(-DL_sponge, 0.0)};

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(emitter_buffer_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
//----------------------------------------------------------------------
//	Define emitter buffer inflow boundary condition
//----------------------------------------------------------------------
class EmitterBufferInflowCondition : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref_;

public:
	EmitterBufferInflowCondition(FluidBody &fluid_body, BodyPartByCell &constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region),
		u_ave_(0), u_ref_(U_f), t_ref_(2.0) {}

	Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
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

	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real run_time_ = GlobalStaticVariables::physical_time_;
		du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);

		return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
	}
};
//----------------------------------------------------------------------
//	Define Observer Particle Generator
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** the measureing particles */
		Vec2d point_coordinate_1(3.0, 5.0);
		Vec2d point_coordinate_2(4.0, 5.0);
		Vec2d point_coordinate_3(5.0, 5.0);

		positions_volumes_.push_back(std::make_pair(point_coordinate_1, 0.0));
		positions_volumes_.push_back(std::make_pair(point_coordinate_2, 0.0));
		positions_volumes_.push_back(std::make_pair(point_coordinate_3, 0.0));
	}
};
