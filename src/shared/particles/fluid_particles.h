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
 * @file 	fluid_particles.h
 * @brief 	This is the derived class of base particle.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef FLUID_PARTICLES_H
#define FLUID_PARTICLES_H

#include "base_particles.hpp"
#include "fluid_particles_variable.h"

#include "particle_generator_lattice.h"

namespace SPH
{
	class Fluid;
	class Oldroyd_B_Fluid;
	class CompressibleFluid;

	/**
	 * @class BaseParticles
	 * @brief newtonian fluid particles.
	 */
	// class BaseParticles : public BaseParticles
	// {
	// public:
	// 	StdLargeVec<Real> p_;				 /**< pressure */
	// 	StdLargeVec<Real> drho_dt_;			 /**< density change rate */
	// 	StdLargeVec<Real> rho_sum_;			 /**< density by particle summation */
	// 	StdLargeVec<int> surface_indicator_; /**< free surface indicator */
	// 	Fluid &fluid_;

	// 	BaseParticles(SPHBody &sph_body, Fluid *fluid);
	// 	virtual ~BaseParticles(){};
	// 	/** Initialize particle variables used in fluid particle. */
	// 	virtual void initializeOtherVariables() override;
	// 	/** Return the ptr of this object. */
	// 	virtual BaseParticles *ThisObjectPtr() override { return this; };
	// };

	// /**
	//  * @class ViscoelasticBaseParticles
	//  * @brief Viscoelastic fluid particles.
	//  */
	// class ViscoelasticBaseParticles : public BaseParticles
	// {
	// public:
	// 	StdLargeVec<Matd> tau_;		/**<  elastic stress */
	// 	StdLargeVec<Matd> dtau_dt_; /**<  change rate of elastic stress */
	// 	Oldroyd_B_Fluid &oldroyd_b_fluid_; 

	// 	ViscoelasticBaseParticles(SPHBody &sph_body, Oldroyd_B_Fluid *oldroyd_b_fluid);
	// 	virtual ~ViscoelasticBaseParticles(){};
	// 	/** Initialize particle variables used in viscoelastic fluid particle. */
	// 	virtual void initializeOtherVariables() override;
	// 	/** Return the ptr of this object. */
	// 	virtual ViscoelasticBaseParticles *ThisObjectPtr() override { return this; };
	// };

	// /**
	//  * @class CompressibleBaseParticles
	//  * @brief Compressible fluid particles.
	//  */
	// class CompressibleBaseParticles : public BaseParticles
	// {
	// public:
	// 	StdLargeVec<Vecd> mom_;				/**< momentum */
	// 	StdLargeVec<Vecd> dmom_dt_; 		/**< change rate of momentum */
	// 	StdLargeVec<Vecd> dmom_dt_prior_;
	// 	StdLargeVec<Real> E_;	  			/**< total energy per unit volume */
	// 	StdLargeVec<Real> dE_dt_; 			/**< change rate of total energy */
	// 	StdLargeVec<Real> dE_dt_prior_;
	// 	CompressibleFluid &compressible_fluid_;

	// 	CompressibleBaseParticles(SPHBody &sph_body, CompressibleFluid *compressible_fluid);
	// 	virtual ~CompressibleBaseParticles(){};
	// 	/** Initialize particle variables used in compressible fluid particle. */
	// 	virtual void initializeOtherVariables() override;
	// 	/** Return the ptr of this object. */
	// 	virtual CompressibleBaseParticles *ThisObjectPtr() override { return this; };
	// };
}
#endif // FLUID_PARTICLES_H