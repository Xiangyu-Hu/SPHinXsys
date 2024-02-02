/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	solid_particles.h
 * @brief 	This is the derived class of base particles.
 * @author	Chi ZHang, Dong Wu and Xiangyu Hu
 */

#ifndef POROUS_SOLID_PARTICLES_H
#define POROUS_SOLID_PARTICLES_H

#include "solid_particles.h"
#include "porous_media_solid.h"
 
namespace SPH
{
namespace multi_species_continuum
{	
 
/**
 * @class PorousMediaParticles
 * @brief A group of particles with elastic body particle data.
 */
	class PorousMediaParticles : public ElasticSolidParticles
	{
	public:
		PorousMediaParticles(SPHBody &body, PorousMediaSolid *porous_solid);
		PorousMediaSolid &porous_solid_;
		virtual ~PorousMediaParticles() {};

		StdLargeVec<Vecd> fluid_velocity_;       /**< fluid velocity */
		StdLargeVec<Vecd> relative_fluid_flux_;  /**<   fluid flux through the boundary of solid unit */
		StdLargeVec<Real> fluid_saturation_;     /**< fluid content quantity in solid */


		StdLargeVec<Real> Vol_update_;             /**< solid volume */
		StdLargeVec<Real> dfluid_mass_dt_;         /**<  fluid mass time gradient in solid */
		StdLargeVec<Real> fluid_mass_;             /**< fluid mass in solid */
		StdLargeVec<Real> total_mass_;             /**< total mass containing fluid mass and solid mass */
		StdLargeVec<Vecd> total_momentum_;         /**< total momentum consists of fluid momentum and solid momentum */
		StdLargeVec<Vecd> dtotal_momentum_dt_;   /**< total momentum time gradient */
		StdLargeVec<Matd> outer_fluid_velocity_relative_fluid_flux_;  /**< outer product of fluid veolcity and fluid flux */
		StdLargeVec<Matd>  Stress_;                            /**< Cauchy stress on solid */
 
		virtual void initializeOtherVariables() override;
		virtual ElasticSolidParticles *ThisObjectPtr() override { return this; };
	};
}
} // namespace SPH
#endif // POROUS_SOLID_PARTICLES_H
