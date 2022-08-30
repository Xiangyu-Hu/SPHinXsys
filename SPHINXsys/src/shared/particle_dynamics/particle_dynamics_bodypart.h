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
 * @file 	particle_dynamics_bodypart.h
 * @brief 	Dynamics for bodypart.
 * The dynamics is constrained to a part of the body,
 * such as in a subregion or on the surface of the body.
 * The particles of a body part can be defined in an Eulerian or Lagrangian fashion.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_BODYPART_H
#define PARTICLE_DYNAMICS_BODYPART_H

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

#include "base_body_part.h"
namespace SPH
{
	/** Body part iterators by particle. sequential computing. */
	void PartIteratorByParticle(const IndexVector &body_part_particles, const ParticleFunctor &particle_functor, Real dt = 0.0);
	/** Body part iterators by particle. parallel computing. */
	void PartIteratorByParticle_parallel(const IndexVector &body_part_particles, const ParticleFunctor &particle_functor, Real dt = 0.0);
	/** Body part iterators by cell. sequential computing. */
	void PartIteratorByCell(const CellLists &body_part_cells, const ParticleFunctor &particle_functor, Real dt = 0.0);
	/** Body part iterators by cell. parallel computing. */
	void PartIteratorByCell_parallel(const CellLists &body_part_cells, const ParticleFunctor &particle_functor, Real dt = 0.0);

	/**
	 * @class PartDynamicsByParticle
	 * @brief Abstract class for imposing body part dynamics by particles.
	 * That is the constrained particles will be the same
	 * during the simulation.
	 */
	class PartDynamicsByParticle : public ParticleDynamics<void>
	{
	public:
		PartDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part);
		virtual ~PartDynamicsByParticle(){};

	protected:
		IndexVector &body_part_particles_;
	};

	/**
	 * @class PartInteractionDynamicsByParticle
	 * @brief Abstract class for particle interaction involving in a body part.
	 */
	class PartInteractionDynamicsByParticle : public PartDynamicsByParticle
	{
	public:
		PartInteractionDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part);
		virtual ~PartInteractionDynamicsByParticle(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		ParticleFunctor functor_interaction_;
		virtual void Interaction(size_t index_i, Real dt = 0.0) = 0;
	};

	/**
	 * @class PartInteractionDynamicsByParticleWithUpdate
	 * @brief Abstract class for particle interaction involving in a body part with an extra update step.
	 */
	class PartInteractionDynamicsByParticleWithUpdate : public PartInteractionDynamicsByParticle
	{
	public:
		PartInteractionDynamicsByParticleWithUpdate(SPHBody &sph_body, BodyPartByParticle &body_part);
		virtual ~PartInteractionDynamicsByParticleWithUpdate(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_update_;
	};

	/**
	 * @class PartInteractionDynamicsByParticleWithUpdate
	 * @brief Abstract class for particle interaction involving in a body part with an extra update step.
	 */
	class PartInteractionDynamicsByParticle1Level : public PartInteractionDynamicsByParticleWithUpdate
	{
	public:
		PartInteractionDynamicsByParticle1Level(SPHBody &sph_body, BodyPartByParticle &body_part);
		virtual ~PartInteractionDynamicsByParticle1Level(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Initialization(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_initialization_;
	};

	/**
	 * @class PartDynamicsByCell
	 * @brief Abstract class for imposing Eulerian constrain to a body.
	 * The constrained particles are in the tagged cells .
	 */
	class PartDynamicsByCell : public ParticleDynamics<void>
	{
	public:
		PartDynamicsByCell(SPHBody &sph_body, BodyPartByCell &body_part);
		virtual ~PartDynamicsByCell(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		CellLists &body_part_cells_;
		ParticleFunctor functor_update_;
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
	};
}
#endif // PARTICLE_DYNAMICS_BODYPART_H