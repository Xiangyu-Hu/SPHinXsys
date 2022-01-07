/**
* @file 	fsi2.h
* @brief 	This is the case file for the test of fluid - structure interaction.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
*/

#ifndef FSI2_CASE_H
#define FSI2_CASE_H

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 11.0;							/**< Channel length. */
Real DH = 4.1;							/**< Channel height. */
Real resolution_ref = 0.1;				/**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;			/**< Boundary width, determined by specific layer of boundary particles. */
Vec2d insert_circle_center(2.0, 2.0);	/**< Location of the cylinder center. */
Real insert_circle_radius = 0.5;		/**< Radius of the cylinder. */
Real bh = 0.4 * insert_circle_radius;	/**< Height of the beam. */
Real bl = 7.0 * insert_circle_radius;	/**< Length of the beam. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;											  /**< Density. */
Real U_f = 1.0;												  /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;										  /**< Speed of sound. */
Real Re = 100.0;											  /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 10.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e3;	/**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, DH));
	water_block_shape.push_back(Vecd(DL, DH));
	water_block_shape.push_back(Vecd(DL, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

	return water_block_shape;
}
/** create a beam shape */
Real hbh = bh / 2.0;
Vec2d BLB(insert_circle_center[0], insert_circle_center[1] - hbh);
Vec2d BLT(insert_circle_center[0], insert_circle_center[1] + hbh);
Vec2d BRB(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] - hbh);
Vec2d BRT(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] + hbh);
std::vector<Vecd> createBeamShape()
{
	std::vector<Vecd> beam_shape;
	beam_shape.push_back(BLB);
	beam_shape.push_back(BLT);
	beam_shape.push_back(BRT);
	beam_shape.push_back(BRB);
	beam_shape.push_back(BLB);

	return beam_shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

	return inner_wall_shape;
}
//----------------------------------------------------------------------
//	Define case dependent bodies material, constraint and boundary conditions.
//----------------------------------------------------------------------
/** Fluid body definition */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(createBeamShape(), ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/* Definition of the solid body. */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/** Definition of the inserted body as a elastic structure. */
class InsertedBody : public SolidBody
{
public:
	InsertedBody(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 2.0))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createBeamShape(), ShapeBooleanOps::add);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};
/** create the beam base for constrain . */
MultiPolygon createBeamBaseShape()
{
	MultiPolygon multi_polygon;
	multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(createBeamShape(), ShapeBooleanOps::sub);
	return multi_polygon;
}
/** create a inflow buffer shape. */
MultiPolygon createInflowBufferShape()
{
	std::vector<Vecd> inflow_buffer_shape;
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, DH));
	inflow_buffer_shape.push_back(Vecd(0.0, DH));
	inflow_buffer_shape.push_back(Vecd(0.0, 0.0));
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
/** Case dependent inflow boundary condition. */
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(FluidBody &fluid_body, BodyPartByCell &constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region),
		  u_ave_(0), u_ref_(1.0), t_ref(2.0) {}
	Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0)
		{
			u = 6.0 * u_ave_ * position[1] * (DH - position[1]) / DH / DH;
			v = 0.0;
		}
		return Vecd(u, v);
	}
	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;
	}
};
/** beam observer particle generator */
class BeamObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	BeamObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(0.5 * (BRT + BRB), 0.0));
	}
};
/** fluid observer particle generator */
class FluidObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	FluidObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** A line of measuring points at the entrance of the channel. */
		size_t number_observation_points = 21;
		Real range_of_measure = DH - resolution_ref * 4.0;
		Real start_of_measure = resolution_ref * 2.0;
		/** the measureing particles */
		for (size_t i = 0; i < number_observation_points; ++i)
		{
			Vec2d point_coordinate(0.0, range_of_measure * (Real)i / (Real)(number_observation_points - 1) + start_of_measure);
			positions_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
		}
	}
};

#endif //FSI2_CASE_H
