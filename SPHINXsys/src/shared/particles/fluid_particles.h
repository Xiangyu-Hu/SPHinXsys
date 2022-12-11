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
	 * @class FluidParticles
	 * @brief newtonian fluid particles.
	 */
	class FluidParticles : public BaseParticles
	{
	public:
		StdLargeVec<Real> p_;				 /**< pressure */
		StdLargeVec<Real> drho_dt_;			 /**< density change rate */
		StdLargeVec<Real> rho_sum_;			 /**< number density */
		StdLargeVec<int> surface_indicator_; /**< free surface indicator */
		Fluid &fluid_;

		FluidParticles(SPHBody &sph_body, Fluid *fluid);
		virtual ~FluidParticles(){};
		/** Initialize particle variables used in fluid particle. */
		virtual void initializeOtherVariables() override;
		/** Return the ptr of this object. */
		virtual FluidParticles *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class ViscoelasticFluidParticles
	 * @brief Viscoelastic fluid particles.
	 */
	class ViscoelasticFluidParticles : public FluidParticles
	{
	public:
		StdLargeVec<Matd> tau_;		/**<  elastic stress */
		StdLargeVec<Matd> dtau_dt_; /**<  change rate of elastic stress */
		Oldroyd_B_Fluid &oldroyd_b_fluid_; 

		ViscoelasticFluidParticles(SPHBody &sph_body, Oldroyd_B_Fluid *oldroyd_b_fluid);
		virtual ~ViscoelasticFluidParticles(){};
		/** Initialize particle variables used in viscoelastic fluid particle. */
		virtual void initializeOtherVariables() override;
		/** Return the ptr of this object. */
		virtual ViscoelasticFluidParticles *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class CompressibleFluidParticles
	 * @brief Compressible fluid particles.
	 */
	class CompressibleFluidParticles : public FluidParticles
	{
	public:
		StdLargeVec<Vecd> mom_;				/**< momentum */
		StdLargeVec<Vecd> dmom_dt_; 		/**< change rate of momentum */
		StdLargeVec<Vecd> dmom_dt_prior_;
		StdLargeVec<Real> E_;	  			/**< total energy per unit volume */
		StdLargeVec<Real> dE_dt_; 			/**< change rate of total energy */
		StdLargeVec<Real> dE_dt_prior_;
		CompressibleFluid &compressible_fluid_;

		CompressibleFluidParticles(SPHBody &sph_body, CompressibleFluid *compressible_fluid);
		virtual ~CompressibleFluidParticles(){};
		/** Initialize particle variables used in compressible fluid particle. */
		virtual void initializeOtherVariables() override;
		/** Return the ptr of this object. */
		virtual CompressibleFluidParticles *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class WeaklyCompressibleFluidParticles
	 * @brief WeaklyCompressible fluid particles.
	 */
	class WeaklyCompressibleFluidParticles : public FluidParticles
	{
	public:
		StdLargeVec<Real> dmass_dt_;	  /**< mass change rate */
		StdLargeVec<Vecd> mom_;			  /**< momentum */
		StdLargeVec<Vecd> dmom_dt_;		  /**< change rate of momentum */
		StdLargeVec<Vecd> dmom_dt_prior_; /**< other, such as gravity and viscous, accelerations, cause momentum loss */

		WeaklyCompressibleFluidParticles(SPHBody &sph_body, Fluid *fluid);
		virtual ~WeaklyCompressibleFluidParticles(){};

		/** Initialize particle variables used in weakly-compressible fluid particle. */
		virtual void initializeOtherVariables() override;
		/** Return the ptr of this object. */
		virtual WeaklyCompressibleFluidParticles *ThisObjectPtr() override { return this; };
	};
}
#endif // FLUID_PARTICLES_H