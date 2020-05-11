/**
* @file 	fsi2_case.h
* @brief 	This is the case file for the test of fliud - structure interaction.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
* @version 0.1
*/

#pragma once

#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters.
 */
Real DL = 11.0; 					/**< Channel length. */
Real DH = 4.1; 						/**< Channel height. */
Real particle_spacing_ref = 0.1; 	/**< Initial reference particle spacing. */
Real DLsponge = particle_spacing_ref * 20.0;	/**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0; 		/**< Boundary width, determined by specific layer of boundary particles. */
Vec2d insert_circle_center(2.0, 2.0);		/**< Location of the cylinder center. */
Real insert_circle_radius = 0.5;			/**< Radius of the cylinder. */
Real bh = 0.4 * insert_circle_radius;			/**< Height of the beam. */
Real bl = 7.0 * insert_circle_radius;			/**< Length of the beam. */
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		/**< Density. */
Real U_f = 1.0;			/**< Cheractristic velocity. */
Real c_f = 10.0 * U_f;	/**< Speed of sound. */
Real Re = 100.0;		/**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re;	/**< Dynamics visocisty. */
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
std::vector<Point> CreatWaterBlockShape()
{
	//geometry
	std::vector<Point> water_block_shape;
	water_block_shape.push_back(Point(-DLsponge, 0.0));
	water_block_shape.push_back(Point(-DLsponge, DH));
	water_block_shape.push_back(Point(DL, DH));
	water_block_shape.push_back(Point(DL, 0.0));
	water_block_shape.push_back(Point(-DLsponge, 0.0));

	return water_block_shape;
}
/** create a water block buffer shape. */
std::vector<Point> CreatInflowBufferShape()
{
	std::vector<Point> inlfow_buffer_shape;
	inlfow_buffer_shape.push_back(Point(-DLsponge, 0.0));
	inlfow_buffer_shape.push_back(Point(-DLsponge, DH));
	inlfow_buffer_shape.push_back(Point(0.0, DH));
	inlfow_buffer_shape.push_back(Point(0.0, 0.0));
	inlfow_buffer_shape.push_back(Point(-DLsponge, 0.0));

	return inlfow_buffer_shape;
}
/** create a beam shape */
Real hbh = bh / 2.0;
Vec2d BLB(insert_circle_center[0], insert_circle_center[1] - hbh);
Vec2d BLT(insert_circle_center[0], insert_circle_center[1] + hbh);
Vec2d BRB(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] - hbh);
Vec2d BRT(insert_circle_center[0] + insert_circle_radius + bl, insert_circle_center[1] + hbh);
std::vector<Point> CreatBeamShape()
{
	std::vector<Point> beam_shape;
	beam_shape.push_back(BLB);
	beam_shape.push_back(BLT);
	beam_shape.push_back(BRT);
	beam_shape.push_back(BRB);
	beam_shape.push_back(BLB);

	return beam_shape;
}
/** create outer wall shape */
std::vector<Point> CreatOuterWallShape()
{
	std::vector<Point> outer_wall_shape;
	outer_wall_shape.push_back(Point(-DLsponge - BW, -BW));
	outer_wall_shape.push_back(Point(-DLsponge - BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, -BW));
	outer_wall_shape.push_back(Point(-DLsponge - BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Point> CreatInnerWallShape()
{
	std::vector<Point> inner_wall_shape;
	inner_wall_shape.push_back(Point(-DLsponge - 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Point(-DLsponge - 2.0 * BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0 * BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Point(-DLsponge - 2.0 * BW, 0.0));

	return inner_wall_shape;
}
/**
 * @brief Define case dependent bodies matertial, constraint and boundary conditions.
 */
 /** Fluid body definition */
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& system, string body_name,	int refinement_level, ParticlesGeneratorOps op)
		: FluidBody(system, body_name, refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> water_bock_shape = CreatWaterBlockShape();
		body_region_.add_geometry(new Geometry(water_bock_shape), RegionBooleanOps::add);
		/** Geomerty definition. */
		body_region_.add_geometry(new Geometry(insert_circle_center, insert_circle_radius, 100), RegionBooleanOps::sub);
		std::vector<Point> beam_shape = CreatBeamShape();
		body_region_.add_geometry(new Geometry(beam_shape), RegionBooleanOps::sub);
		/** Finalize the geometry definition and correspoding opertation. */
		body_region_.done_modeling();
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
		/** supplementrary material paramters derived from basic parameters. */
		assignDerivedMaterialParameters();
	}
};
/* Definition of the solid body. */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem& system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> outer_wall_shape = CreatOuterWallShape();
		std::vector<Point> inner_wall_shape = CreatInnerWallShape();
		body_region_.add_geometry(new Geometry(outer_wall_shape), RegionBooleanOps::add);
		body_region_.add_geometry(new Geometry(inner_wall_shape), RegionBooleanOps::sub);
		/** Finalize the geometry definition and correspoding opertation. */
		body_region_.done_modeling();
	}
};
/** Definition of the inserted body as a elastic structure. */
class InsertedBody : public SolidBody
{
public:
	InsertedBody(SPHSystem& system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		/** Geomerty definition. */
		std::vector<Point> beam_shape = CreatBeamShape();
		Geometry* circle_geometry = new Geometry(insert_circle_center, insert_circle_radius, 100);
		body_region_.add_geometry(circle_geometry, RegionBooleanOps::add);
		Geometry* beam_geometry = new Geometry(beam_shape);
		body_region_.add_geometry(beam_geometry, RegionBooleanOps::add);
		/** Finalize the geometry definition and correspoding opertation. */
		body_region_.done_modeling();
	}
};
/** the material for insert body. */
class InsertBodyMaterial : public LinearElasticSolid
{
public:
	InsertBodyMaterial() : LinearElasticSolid()
	{
		rho_0_ = rho0_s;
		E_0_ = Youngs_modulus;
		nu_ = poisson;

		assignDerivedMaterialParameters();
	}
};
/** constraint the cylinder part of the insert body. */
class BeamBase : public BodyPartByParticle
{
public:
	BeamBase(SolidBody* solid_body, string constrianed_region_name)
		: BodyPartByParticle(solid_body, constrianed_region_name)
	{
		/** Geomerty definition. */
		std::vector<Point> beam_shape = CreatBeamShape();
		Geometry* circle_geometry = new Geometry(insert_circle_center, insert_circle_radius, 100);
		body_part_region_.add_geometry(circle_geometry, RegionBooleanOps::add);
		Geometry* beam_gemetry = new Geometry(beam_shape);
		body_part_region_.add_geometry(beam_gemetry, RegionBooleanOps::sub);
		body_part_region_.done_modeling();
		/**  Tag the constrained particle. */
		TagBodyPartParticles();
	}
};
/** inflow buffer */
class InflowBuffer : public BodyPartByCell
{
public:
	InflowBuffer(FluidBody* fluid_body, string constrianed_region_name)
		: BodyPartByCell(fluid_body, constrianed_region_name)
	{
		/** Geomerty definition. */
		std::vector<Point> inflow_buffer_shape = CreatInflowBufferShape();
		body_part_region_.add_geometry(new Geometry(inflow_buffer_shape), RegionBooleanOps::add);
		/** Finalize the geometry definition and correspoding opertation. */
		body_part_region_.done_modeling();
		//tag the constrained particle
		TagBodyPartCells();
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
	Vecd GetInflowVelocity(Vecd& position, Vecd& velocity)
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0) {
			u = 6.0 * u_ave_ * position[1] * (DH - position[1]) / DH / DH;
			v = 0.0;
		}
		return Vecd(u, v);
	}
	void PrepareConstraint() override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(pi * run_time / t_ref)) : u_ref_;
	}
};
/** fluid observer body */
class BeamObserver : public FictitiousBody
{
public:
	BeamObserver(SPHSystem& system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{
		/** the measuring particle with zero volume */
		body_input_points_volumes_.push_back(make_pair(0.5 * (BRT + BRB), 0.0));
	}
};
/** an observer body to measure the flow profile */
class FluidObserver : public FictitiousBody
{
public:
	FluidObserver(SPHSystem& system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: FictitiousBody(system, body_name, refinement_level, 1.3, op)
	{
		/** A line of measuring points at the entrance of the channel. */
		size_t number_observation_pionts = 21;
		Real range_of_measure = DH - particle_spacing_ref * 4.0;
		Real start_of_measure = particle_spacing_ref * 2.0;
		/** the measureing particles */
		for (size_t i = 0; i < number_observation_pionts; ++i) {
			Vec2d point_coordinate(0.0, range_of_measure * Real(i) / Real(number_observation_pionts - 1) + start_of_measure);
			body_input_points_volumes_.push_back(make_pair(point_coordinate, 0.0));
		}
	}
};