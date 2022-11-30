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
* @brief 	Here, we test the fluid-shell interaction model with dam break impact thin structure. 
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
Real DL = 0.8; 								     	/**< Tank length. */
Real DH = 0.6; 							         	/**< Tank height. */
Real L_W = 0.2;                                   	/**< water width. */
Real L_H = 0.3;                                    	/**< water depth. */
Real Gate_x = 0.6;						  	 		/**< X-axis of gate */
Real Gate_width = 0.004;							/**< Width of the gate. */
Real Gate_height = 0.1;						 		/**< Height of the gate. */
Real particle_spacing_ref = Gate_width; 		/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4.0; 		   		/**< Extending width for BCs. */
Real thickness = Gate_width;
int particle_number_mid_surface = int((Gate_height) / particle_spacing_ref);
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 997.0;						/**< Reference density of fluid. */
Real gravity_g = 9.8; 					    /**< Value of gravity. */
Real U_max = 2.0 * sqrt(gravity_g * L_H);	/**< Characteristic velocity. */
Real c_f = 10.0 * U_max;					/**< Reference sound speed. */
Real mu_f = 8.9e-4;							/**< Dynamics viscosity. */
/**
 * @brief Material properties of the elastic gate.
 */
Real rho0_s = 1161.54;				/**< Reference density of gate. */
Real Youngs_modulus = 3.5e6;		/**< reference Youngs modulus. */
Real poisson = 0.45; 				/**< Poisson ratio. */
Real physical_viscosity = 1000.0; 	/**< physical damping, here we choose the same value as numerical viscosity. */
/** create a water block shape */
StdVec<Vecd> water_block_shape{
	Vecd(0.0, 0.0), Vecd(0.0, L_H), Vecd(L_W, L_H), Vecd(L_W, 0.0), Vecd(0.0, 0.0)};
/** Inner wall shape */
StdVec<Vecd> inner_wall_shape{
    Vecd(0.0, 0.0),
	Vecd(0.0, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(0.0, 0.0)
};
/** Outer wall shape */
StdVec<Vecd> outer_wall_shape{
	Vecd(-BW, -BW),
	Vecd(-BW, DH + BW),
	Vecd(DL + BW, DH + BW),
	Vecd(DL + BW, -BW),
	Vecd(-BW, -BW)
};
/* Case-dependent geometries.*/
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};
/* Case-dependent geometries.*/
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
/** Particle generator and constraint boundary for shell baffle. */
class ShellBaffleParticleGenerator : public SurfaceParticleGenerator
{
public:
	explicit ShellBaffleParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
	virtual void initializeGeometricVariables() override
	{
		for (int i = 0; i < particle_number_mid_surface; i++)
		{
			Real x = Gate_x - 0.5 * particle_spacing_ref;
			Real y = -0.5 * particle_spacing_ref + Real(i) * particle_spacing_ref;
            initializePositionAndVolumetricMeasure(Vecd(x, y), particle_spacing_ref);
            Vec2d normal_direction = Vec2d(1.0, 0);
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
		if (base_particles_.pos_[index_i][1] < 0.0)
		{
			body_part_particles_.push_back(index_i);
		}
	};
};

/**
 * @brief Define the observer body.
 */
StdVec<Vecd> baffle_disp_probe_location{Vecd(Gate_x - 0.5 * particle_spacing_ref, 0.0785), Vecd(Gate_x - 0.5 * particle_spacing_ref, 0.065), Vecd(Gate_x - 0.5 * particle_spacing_ref, 0.04)};