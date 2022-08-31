/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
* @file 	particle_dynamics_algorithms.h
* @brief 	This is the classes for algorithms particle dynamics.
* @detail	Generally, there are four types dynamics. One is without particle interaction.
*			One is with particle interaction within a body. One is with particle interaction
*			between a center body and other contacted bodies.
* 			Still another is the combination of the last two.
*			For the first dynamics, there is also reduce dynamics
*			which carries reduced operations through the particles of the body.

* @author	Chi ZHang and Xiangyu Hu
*/

#ifndef PARTICLE_DYNAMICS_ALGORITHMS_H
#define PARTICLE_DYNAMICS_ALGORITHMS_H

#include "particle_iterators.h"
#include "base_local_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH
{
	/**
	 * @class InteractionDynamics
	 * @brief This is the class for particle interaction with other particles
	 */
	class InteractionDynamics : public ParticleDynamics<void>
	{
	public:
		explicit InteractionDynamics(SPHBody &sph_body)
			: ParticleDynamics<void>(sph_body),
			  functor_interaction_(std::bind(&InteractionDynamics::Interaction,
											 this, _1, _2)){};
		virtual ~InteractionDynamics(){};

		/** pre process such as update ghost state */
		StdVec<BaseDynamics<void> *> pre_processes_;
		/** post process such as impose constraint */
		StdVec<BaseDynamics<void> *> post_processes_;

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		friend class CombinedInteractionDynamics;
		virtual void Interaction(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_interaction_;
	};

	/**
	 * @class CombinedInteractionDynamics
	 * @brief This is the class for combining several interactions dynamics,
	 * which share the particle loop but are independent from each other,
	 * aiming to increase computing intensity under the data caching environment
	 */
	class CombinedInteractionDynamics : public InteractionDynamics
	{
	public:
		explicit CombinedInteractionDynamics(InteractionDynamics &dynamics_a, InteractionDynamics &dynamics_b);
		virtual ~CombinedInteractionDynamics(){};

	protected:
		InteractionDynamics &dynamics_a_, &dynamics_b_;
		virtual void setupDynamics(Real dt = 0.0) override;
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class InteractionDynamicsWithUpdate
	 * @brief This class includes an interaction and a update steps
	 */
	class InteractionDynamicsWithUpdate : public InteractionDynamics
	{
	public:
		explicit InteractionDynamicsWithUpdate(SPHBody &sph_body)
			: InteractionDynamics(sph_body),
			  functor_update_(std::bind(&InteractionDynamicsWithUpdate::Update,
										this, _1, _2)) {}
		virtual ~InteractionDynamicsWithUpdate(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_update_;
	};

	/**
	 * @class ParticleDynamics1Level
	 * @brief This class includes an initialization, an interaction and a update steps
	 */
	class ParticleDynamics1Level : public InteractionDynamicsWithUpdate
	{
	public:
		explicit ParticleDynamics1Level(SPHBody &sph_body)
			: InteractionDynamicsWithUpdate(sph_body),
			  functor_initialization_(std::bind(&ParticleDynamics1Level::Initialization,
												this, _1, _2)) {}
		virtual ~ParticleDynamics1Level(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Initialization(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_initialization_;
	};

	/**
	 * @class InteractionDynamicsSplitting
	 * @brief This is for the splitting algorithm
	 */
	class InteractionDynamicsSplitting : public InteractionDynamics
	{
	public:
		explicit InteractionDynamicsSplitting(SPHBody &sph_body);
		virtual ~InteractionDynamicsSplitting(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		RealBody &real_body_;
		SplitCellLists &split_cell_lists_;
	};

	//----------------------------------------------------------------------
	//		New version particle dynamics base classes.
	//		Aiming to use template on local dynamics so that
	//		it can be used in different dynamics.
	//----------------------------------------------------------------------

	/**
	 * @class SimpleDynamics
	 * @brief Simple particle dynamics without considering particle interaction
	 */
	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class SimpleDynamics : public LocalDynamicsType, public BaseDynamics<void>
	{
		DynamicsRange &dynamics_range_;

	public:
		template <typename... Args>
		SimpleDynamics(DynamicsRange &dynamics_range, Args &&...args)
			: LocalDynamicsType(dynamics_range, std::forward<Args>(args)...),
			  BaseDynamics<void>(), dynamics_range_(dynamics_range){};
		virtual ~SimpleDynamics(){};

		virtual void exec(Real dt = 0.0) override
		{
			this->setBodyUpdated();
			this->setupDynamics(dt);
			particle_for(
				dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->update(i, delta); },
				dt);
		};

		virtual void parallel_exec(Real dt = 0.0) override
		{
			this->setBodyUpdated();
			this->setupDynamics(dt);
			particle_parallel_for(
				dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->update(i, delta); },
				dt);
		};
	};

	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class ReduceDynamics : public LocalDynamicsType,
						   public BaseDynamics<typename LocalDynamicsType::ReduceReturnType>

	{
		using ReturnType = typename LocalDynamicsType::ReduceReturnType;

	protected:
		DynamicsRange &dynamics_range_;

	public:
		template <typename... Args>
		ReduceDynamics(DynamicsRange &dynamics_range, Args &&...args)
			: LocalDynamicsType(dynamics_range, std::forward<Args>(args)...),
			  BaseDynamics<ReturnType>(), dynamics_range_(dynamics_range){};
		virtual ~ReduceDynamics(){};

		using ReduceReturnType = ReturnType;
		std::string QuantityName() { return this->quantity_name_; };
		std::string DynamicsRangeName() { return dynamics_range_.getName(); };

		virtual ReturnType exec(Real dt = 0.0) override
		{
			this->setupDynamics(dt);
			ReturnType temp = particle_reduce(
				dynamics_range_.LoopRange(), this->Reference(), this->getOperation(),
				[&](size_t i, Real delta) -> ReturnType
				{ return this->reduce(i, delta); },
				dt);
			return this->outputResult(temp);
		};

		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			this->setupDynamics(dt);
			ReturnType temp = particle_parallel_reduce(
				dynamics_range_.LoopRange(), this->Reference(), this->getOperation(),
				[&](size_t i, Real delta) -> ReturnType
				{ return this->reduce(i, delta); },
				dt);
			return this->outputResult(temp);
		};
	};

	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class ReduceDynamicsAverage : public ReduceDynamics<LocalDynamicsType, DynamicsRange>
	{
		using ReturnType = typename LocalDynamicsType::ReduceReturnType;
		ReturnType outputAverage(ReturnType sum, size_t size_of_loop_range)
		{
			return sum / Real(size_of_loop_range);
		}

	public:
		template <typename... Args>
		ReduceDynamicsAverage(DynamicsRange &dynamics_range, Args &&...args)
			: ReduceDynamics<LocalDynamicsType, DynamicsRange>(dynamics_range, std::forward<Args>(args)...){};
		virtual ~ReduceDynamicsAverage(){};

		virtual ReturnType exec(Real dt = 0.0) override
		{
			ReturnType sum = ReduceDynamics<LocalDynamicsType, DynamicsRange>::exec(dt);
			return outputAverage(sum, this->dynamics_range_.SizeOfLoopRange());
		};

		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			ReturnType sum = ReduceDynamics<LocalDynamicsType, DynamicsRange>::parallel_exec(dt);
			return outputAverage(sum, this->dynamics_range_.SizeOfLoopRange());
		};
	};
}
#endif // PARTICLE_DYNAMICS_ALGORITHMS_H