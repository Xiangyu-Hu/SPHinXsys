/**
* @file 	fsi2_case.h
* @brief 	This is the case file for the test of fluid - structure interaction.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
*/

#ifndef FSI2_CASE_H
#define FSI2_CASE_H

#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters.
 */
Real DL = 11.0; 								/**< Channel length. */
Real DH = 4.1; 									/**< Channel height. */
Real resolution_ref = 0.1; 						/**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0;	/**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0; 			/**< Boundary width, determined by specific layer of boundary particles. */
Vec2d insert_circle_center(2.0, 2.0);			/**< Location of the cylinder center. */
Real insert_circle_radius = 0.5;				/**< Radius of the cylinder. */
Real bh = 0.4 * insert_circle_radius;			/**< Height of the beam. */
Real bl = 7.0 * insert_circle_radius;			/**< Length of the beam. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		/**< Density. */
Real U_f = 1.0;			/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;	/**< Speed of sound. */
Real Re = 100.0;		/**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re;	/**< Dynamics viscosity. */
/**
 * @brief Material properties of the solid,
 */
Real rho0_s = 10.0; 		/**< Reference density.*/
Real poisson = 0.4; 		/**< Poisson ratio.*/
Real Ae = 1.4e3; 			/**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
/**
* @brief define geometry of SPH bodies
*/
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
/** create a water block buffer shape. */
std::vector<Vecd> CreatInflowBufferShape()
{
	std::vector<Vecd> inlfow_buffer_shape;
	inlfow_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
	inlfow_buffer_shape.push_back(Vecd(-DL_sponge, DH));
	inlfow_buffer_shape.push_back(Vecd(0.0, DH));
	inlfow_buffer_shape.push_back(Vecd(0.0, 0.0));
	inlfow_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

	return inlfow_buffer_shape;
}
/** create a beam shape */
Real hbh = bh / 2.0;
Vec2d BLB(insert_circle_center[0], insert_circle_center[1] - hbh);
Vec2d BLT(insert_circle_center[0], insert_circle_center[1] + hbh);
Vec2d BRB(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] - hbh);
Vec2d BRT(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] + hbh);
std::vector<Vecd> CreatBeamShape()
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
		/** Geomtry definition. */
		std::vector<Vecd> beam_shape = CreatBeamShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);
		body_shape_->addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
		body_shape_->addAPolygon(beam_shape, ShapeBooleanOps::sub);
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
/* Definition of the solid body. */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem& system, std::string body_name)
		: SolidBody(system, body_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> outer_wall_shape = createOuterWallShape();
		std::vector<Vecd> inner_wall_shape = createInnerWallShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
/** Definition of the inserted body as a elastic structure. */
class InsertedBody : public SolidBody
{
public:
	InsertedBody(SPHSystem& system, std::string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(1.15, 2.0))
	{
		/** Geomtry definition. */
		std::vector<Vecd> beam_shape = CreatBeamShape();
		ComplexShape original_body_shape;
		original_body_shape.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
		original_body_shape.addAPolygon(beam_shape, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);

	}
};
/** the material for insert body. */
class InsertBodyMaterial : public LinearElasticSolid
{
public:
	InsertBodyMaterial() : LinearElasticSolid()
	{
		rho0_ = rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};
/** constraint the cylinder part of the insert body. */
class BeamBase : public BodyPartByParticle
{
public:
	BeamBase(SolidBody* solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> beam_shape = CreatBeamShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
		body_part_shape_->addAPolygon(beam_shape, ShapeBooleanOps::sub);

		/**  Tag the constrained particle. */
		tagBodyPart();
	}
};
/** inflow buffer */
class InflowBuffer : public BodyPartByCell
{
public:
	InflowBuffer(FluidBody* fluid_body, std::string constrained_region_name)
		: BodyPartByCell(fluid_body, constrained_region_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> inflow_buffer_shape = CreatInflowBufferShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);

		//tag the constrained particle
		tagBodyPart();
	}
};
/** Case dependent inflow boundary condition. */
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, u_ref_, t_ref;
public:
	ParabolicInflow(FluidBody* fluid_body,
		BodyPartByCell* constrained_region)
		: InflowBoundaryCondition(fluid_body, constrained_region)
	{
		u_ave_ = 0.0;
		u_ref_ = 1.0;
		t_ref = 2.0;
	}
	Vecd getTargetVelocity(Vecd& position, Vecd& velocity)
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0) {
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
/** fluid observer body */
class BeamObserver : public FictitiousBody
{
public:
	BeamObserver(SPHSystem& system, std::string body_name)
		: FictitiousBody(system, body_name, new ParticleAdaptation(1.15, 2.0))
	{
		/** the measuring particle with zero volume */
		body_input_points_volumes_.push_back(std::make_pair(0.5 * (BRT + BRB), 0.0));
	}
};
/** an observer body to measure the flow profile */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem& system, std::string body_name)
		: FictitiousBody(system, body_name)
	{
		/** A line of measuring points at the entrance of the channel. */
		size_t number_observation_pionts = 21;
		Real range_of_measure = DH - resolution_ref * 4.0;
		Real start_of_measure = resolution_ref * 2.0;
		/** the measureing particles */
		for (size_t i = 0; i < number_observation_pionts; ++i) {
			Vec2d point_coordinate(0.0, range_of_measure * Real(i) / Real(number_observation_pionts - 1) + start_of_measure);
			body_input_points_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
		}
	}
};
#endif //FSI2_CASE_H