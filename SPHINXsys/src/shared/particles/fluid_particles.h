/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	fluid_particles.h
 * @brief 	This is the derived class of base particle.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef FLUID_PARTICLES_H
#define FLUID_PARTICLES_H

#include "base_particles.h"
#include "base_particles.hpp"

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
		explicit FluidParticles(SPHBody &sph_body,
								SharedPtr<Fluid> shared_fluid_ptr,
								SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>());
		virtual ~FluidParticles(){};

		StdLargeVec<Real> p_;				 /**< pressure */
		StdLargeVec<Real> drho_dt_;			 /**< density change rate */
		StdLargeVec<Real> rho_sum_;			 /**< number density */
		StdLargeVec<int> surface_indicator_; /**< free surface indicator */

		virtual FluidParticles *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class ViscoelasticFluidParticles
	 * @brief Viscoelastic fluid particles.
	 */
	class ViscoelasticFluidParticles : public FluidParticles
	{
	public:
		explicit ViscoelasticFluidParticles(SPHBody &sph_body,
											SharedPtr<Oldroyd_B_Fluid> shared_oldroyd_b_fluid_ptr,
											SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>());
		virtual ~ViscoelasticFluidParticles(){};

		StdLargeVec<Matd> tau_;		/**<  elastic stress */
		StdLargeVec<Matd> dtau_dt_; /**<  change rate of elastic stress */

		virtual ViscoelasticFluidParticles *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class CompressibleFluidParticles
	 * @brief Compressible fluid particles.
	 */
	class CompressibleFluidParticles : public FluidParticles
	{
	public:
		explicit CompressibleFluidParticles(SPHBody &sph_body,
											SharedPtr<CompressibleFluid> compressiblefluid,
											SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>());
		virtual ~CompressibleFluidParticles(){};

		StdLargeVec<Vecd> mom_;		/**< momentum */
		StdLargeVec<Vecd> dmom_dt_; /**< change rate of momentum */
		StdLargeVec<Vecd> dmom_dt_prior_;
		StdLargeVec<Real> E_;	  /**< total energy per unit volume */
		StdLargeVec<Real> dE_dt_; /**< change rate of total energy */
		StdLargeVec<Real> dE_dt_prior_;

		virtual CompressibleFluidParticles *ThisObjectPtr() override { return this; };
	};

	/**
	 * @class WeaklyCompressibleFluidParticles
	 * @brief WeaklyCompressible fluid particles.
	 */
	class WeaklyCompressibleFluidParticles : public FluidParticles
	{
	public:
		explicit WeaklyCompressibleFluidParticles(SPHBody &sph_body,
											      SharedPtr<Fluid> shared_fluid_ptr,
											      SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>());
		virtual ~WeaklyCompressibleFluidParticles() {};

		StdLargeVec<Real> dmass_dt_;		/**< mass change rate */
		StdLargeVec<Vecd> mom_;             /**< momentum */
		StdLargeVec<Vecd> dmom_dt_;         /**< change rate of momentum */
		StdLargeVec<Vecd> dmom_dt_prior_;	/**< other, such as gravity and viscous, accelerations, cause momentum loss */

		virtual WeaklyCompressibleFluidParticles* ThisObjectPtr() override { return this; };
	};
}
#endif //FLUID_PARTICLES_H