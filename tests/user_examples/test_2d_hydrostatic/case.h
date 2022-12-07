/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
* @file 	case.h
* @brief 	Here, we test the fluid-shell interaction model with hydrostatic fsi. 
* @author 	Chi Zhang
*/

#include "sphinxsys.h"
#define PI 3.1415926
 /**
* @brief Namespace cite here.
*/
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 1.0;							/**< Tank length. */
Real DH = 2.1;							/**< Tank height. */
Real Dam_L = 1.0;						/**< Water block width. */
Real Dam_H = 2.0;						/**< Water block height. */
Real Gate_width = 0.05;					/**< Width of the gate. */
Real particle_spacing_ref = Gate_width; /**< Initial reference particle spacing. 8, 10, 12 */
Real BW = 4.0 * particle_spacing_ref;		  	/**< Extending width for BCs. */
Real thickness = Gate_width;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;  						/**< Reference density of fluid. */
Real gravity_g = 9.81; 						/**< Value of gravity. */
Real U_f = 2.0 * sqrt(Dam_H * gravity_g);	/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;			  			/**< Reference sound speed. */
Real mu_f = rho0_f * U_f * DL / 0.1; 		/**< For damping. */
//----------------------------------------------------------------------
//	Material properties of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 2700.0;				/**< Reference density of gate. */
Real Youngs_modulus = 6.75e10;		/**< reference Youngs modulus. */
Real poisson = 0.49; 				/**< Poisson ratio. */
Real physical_viscosity = 5000.0; 	/**< physical damping, here we choose the same value as numerical viscosity. */
/** create a water block shape */
//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
	std::vector<Vecd> block_shape;
	block_shape.push_back(Vecd(0,0));
	block_shape.push_back(Vecd(0,Dam_H));
	block_shape.push_back(Vecd(Dam_L,Dam_H));
	block_shape.push_back(Vecd(Dam_L,0));
	block_shape.push_back(Vecd(0,0));

	return block_shape;
}
//----------------------------------------------------------------------
//	wall body definition.
//----------------------------------------------------------------------
std::vector<Vecd> createLeftWallShape()
{
	std::vector<Vecd> wall_shape;
	wall_shape.push_back(Vecd(-BW, -BW));
	wall_shape.push_back(Vecd(-BW, DH));
	wall_shape.push_back(Vecd(0.0, DH));
	wall_shape.push_back(Vecd(0.0, -BW));
	wall_shape.push_back(Vecd(-BW, -BW));

	return wall_shape;
}
std::vector<Vecd> createRightWallShape()
{
	std::vector<Vecd> wall_shape;
	wall_shape.push_back(Vecd(DL, -BW));
	wall_shape.push_back(Vecd(DL, DH));
	wall_shape.push_back(Vecd(DL + BW, DH));
	wall_shape.push_back(Vecd(DL + BW, -BW));
	wall_shape.push_back(Vecd(DL, -BW));

	return wall_shape;
}
//----------------------------------------------------------------------
//	Define gate
//----------------------------------------------------------------------
std::vector<Vecd> createFlexibleShape()
{
	std::vector<Vecd> gate_shape;
	gate_shape.push_back(Vec2d(-BW, -Gate_width));
	gate_shape.push_back(Vec2d(-BW, 0.0));
	gate_shape.push_back(Vec2d(Dam_L + BW, 0.0));
	gate_shape.push_back(Vec2d(Dam_L + BW, -Gate_width));
	gate_shape.push_back(Vec2d(-BW, -Gate_width));

	return gate_shape;
}
/* Case-dependent geometries.*/
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
	}
};
/* Case-dependent geometries.*/
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createLeftWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createRightWallShape(), ShapeBooleanOps::add);
	}
};
/** Particle generator and constraint boundary for shell baffle. */
int particle_number_mid_surface = int((0.0 * BW + Dam_L) / particle_spacing_ref) + 2;
class ShellBaffleParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit ShellBaffleParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number_mid_surface; i++)
		{
			Real x = -particle_spacing_ref + (0.5 + Real(i)) * particle_spacing_ref;
			Real y = - 0.5 * thickness;
            initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_ref);
            Vec2d normal_direction = Vec2d(0, 1.0);
            initializeSurfaceProperties(normal_direction, thickness);
		}
	}
};
/** Define the boundary geometry. */
class BoundaryGeometry : public BodyPartByParticle
{
public:
	BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
		: BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	};
	virtual ~BoundaryGeometry(){};

private:
	void tagManually(size_t index_i)
	{
		if (base_particles_.pos_[index_i][0] <= 0.0 || base_particles_.pos_[index_i][0] >= DL)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/**
 * @brief Define the observer body.
 */
StdVec<Vecd> baffle_disp_probe_location{Vecd(0.5 * DL, - 0.5 * Gate_width)};