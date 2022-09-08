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
	 * @class OldInteractionDynamics
	 * @brief This is the class for particle interaction with other particles
	 */
	class OldInteractionDynamics : public ParticleDynamics<void>
	{
	public:
		explicit OldInteractionDynamics(SPHBody &sph_body)
			: ParticleDynamics<void>(sph_body),
			  functor_interaction_(std::bind(&OldInteractionDynamics::Interaction,
											 this, _1, _2)){};
		virtual ~OldInteractionDynamics(){};

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
	 * @class OldInteractionDynamicsWithUpdate
	 * @brief This class includes an interaction and a update steps
	 */
	class OldInteractionDynamicsWithUpdate : public OldInteractionDynamics
	{
	public:
		explicit OldInteractionDynamicsWithUpdate(SPHBody &sph_body)
			: OldInteractionDynamics(sph_body),
			  functor_update_(std::bind(&OldInteractionDynamicsWithUpdate::Update,
										this, _1, _2)) {}
		virtual ~OldInteractionDynamicsWithUpdate(){};

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
	class ParticleDynamics1Level : public OldInteractionDynamicsWithUpdate
	{
	public:
		explicit ParticleDynamics1Level(SPHBody &sph_body)
			: OldInteractionDynamicsWithUpdate(sph_body),
			  functor_initialization_(std::bind(&ParticleDynamics1Level::Initialization,
												this, _1, _2)) {}
		virtual ~ParticleDynamics1Level(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Initialization(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_initialization_;
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
		template <class DerivedDynamicsRange, typename... Args>
		SimpleDynamics(DerivedDynamicsRange &derived_dynamics_range, Args &&...args)
			: LocalDynamicsType(derived_dynamics_range, std::forward<Args>(args)...),
			  BaseDynamics<void>(), dynamics_range_(derived_dynamics_range){};
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

	/**
	 * @class ReduceDynamics
	 * @brief Template class for particle-wise reduce operation
	 */
	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class ReduceDynamics : public LocalDynamicsType,
						   public BaseDynamics<typename LocalDynamicsType::ReduceReturnType>

	{
		using ReturnType = typename LocalDynamicsType::ReduceReturnType;

	protected:
		DynamicsRange &dynamics_range_;

	public:
		template <class DerivedDynamicsRange, typename... Args>
		ReduceDynamics(DerivedDynamicsRange &derived_dynamics_range, Args &&...args)
			: LocalDynamicsType(derived_dynamics_range, std::forward<Args>(args)...),
			  BaseDynamics<ReturnType>(), dynamics_range_(derived_dynamics_range){};
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

	/**
	 * @class ReduceDynamicsAverage
	 * @brief Template class for computing particle-wise averages
	 */
	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class ReduceDynamicsAverage : public ReduceDynamics<LocalDynamicsType, DynamicsRange>
	{
		using ReturnType = typename LocalDynamicsType::ReduceReturnType;
		ReturnType outputAverage(ReturnType sum, size_t size_of_loop_range)
		{
			return sum / Real(size_of_loop_range);
		}

	public:
		template <class DerivedDynamicsRange, typename... Args>
		ReduceDynamicsAverage(DerivedDynamicsRange &derived_dynamics_range, Args &&...args)
			: ReduceDynamics<LocalDynamicsType, DynamicsRange>(derived_dynamics_range, std::forward<Args>(args)...){};
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

	/**
	 * @class BaseInteractionDynamics
	 * @brief This is the base class for particle interaction with other particles
	 */
	template <class LocalDynamicsType>
	class BaseInteractionDynamics : public LocalDynamicsType, public BaseDynamics<void>
	{
	public:
		template <class BodyRelationType, typename... Args>
		BaseInteractionDynamics(BodyRelationType &body_relation, Args &&...args)
			: LocalDynamicsType(body_relation, std::forward<Args>(args)...),
			  BaseDynamics<void>(){};
		virtual ~BaseInteractionDynamics(){};

		/** pre process such as update ghost state */
		StdVec<BaseDynamics<void> *> pre_processes_;
		/** post process such as impose constraint */
		StdVec<BaseDynamics<void> *> post_processes_;

		virtual void runInteractionStep(Real dt) = 0;
		virtual void parallel_runInteractionStep(Real dt) = 0;

		virtual void exec(Real dt = 0.0) override
		{
			this->setBodyUpdated();
			this->setupDynamics(dt);
			runInteractionStep(dt);
		};

		virtual void parallel_exec(Real dt = 0.0) override
		{
			this->setBodyUpdated();
			this->setupDynamics(dt);
			parallel_runInteractionStep(dt);
		};
	};

	/**
	 * @class NewInteractionDynamicsSplit
	 * @brief This is for the splitting algorithm
	 */
	template <class LocalDynamicsType>
	class NewInteractionDynamicsSplit : public BaseInteractionDynamics<LocalDynamicsType>
	{
	protected:
		RealBody &real_body_;
		SplitCellLists &split_cell_lists_;

	public:
		template <class BodyRelationType, typename... Args>
		NewInteractionDynamicsSplit(BodyRelationType &body_relation, Args &&...args)
			: BaseInteractionDynamics<LocalDynamicsType>(body_relation, std::forward<Args>(args)...),
			  real_body_(DynamicCast<RealBody>(this, this->sph_body_)),
			  split_cell_lists_(real_body_.getSplitCellLists())
		{
			real_body_.setUseSplitCellLists();
		};
		virtual ~NewInteractionDynamicsSplit(){};

		virtual void runInteractionStep(Real dt) override
		{
			for (size_t k = 0; k < this->pre_processes_.size(); ++k)
				this->pre_processes_[k]->exec(dt);

			particle_for_split(
				split_cell_lists_,
				[&](size_t i, Real delta)
				{ this->interaction(i, delta); },
				dt);

			for (size_t k = 0; k < this->post_processes_.size(); ++k)
				this->post_processes_[k]->exec(dt);
		}

		virtual void parallel_runInteractionStep(Real dt) override
		{
			for (size_t k = 0; k < this->pre_processes_.size(); ++k)
				this->pre_processes_[k]->parallel_exec(dt);

			particle_parallel_for_split(
				split_cell_lists_,
				[&](size_t i, Real delta)
				{ this->interaction(i, delta); },
				dt);

			for (size_t k = 0; k < this->post_processes_.size(); ++k)
				this->post_processes_[k]->parallel_exec(dt);
		}
	};

	/**
	 * @class InteractionDynamics
	 * @brief This is the class with a single step of particle interaction with other particles
	 */
	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class InteractionDynamics : public BaseInteractionDynamics<LocalDynamicsType>
	{
	protected:
		DynamicsRange &dynamics_range_;

	public:
		template <class BodyRelationType, typename... Args>
		InteractionDynamics(BodyRelationType &body_relation, Args &&...args)
			: BaseInteractionDynamics<LocalDynamicsType>(body_relation, std::forward<Args>(args)...),
			  dynamics_range_(body_relation.getDynamicsRange()){};
		virtual ~InteractionDynamics(){};

		virtual void runInteractionStep(Real dt) override
		{
			for (size_t k = 0; k < this->pre_processes_.size(); ++k)
				this->pre_processes_[k]->exec(dt);

			particle_for(
				dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->interaction(i, delta); },
				dt);

			for (size_t k = 0; k < this->post_processes_.size(); ++k)
				this->post_processes_[k]->exec(dt);
		}

		virtual void parallel_runInteractionStep(Real dt) override
		{
			for (size_t k = 0; k < this->pre_processes_.size(); ++k)
				this->pre_processes_[k]->parallel_exec(dt);

			particle_parallel_for(
				dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->interaction(i, delta); },
				dt);

			for (size_t k = 0; k < this->post_processes_.size(); ++k)
				this->post_processes_[k]->parallel_exec(dt);
		}
	};

	/**
	 * @class InteractionDynamicsWithUpdate
	 * @brief This class includes an interaction and a update steps
	 */
	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class InteractionDynamicsWithUpdate : public InteractionDynamics<LocalDynamicsType, DynamicsRange>
	{
	public:
		template <class BodyRelationType, typename... Args>
		InteractionDynamicsWithUpdate(BodyRelationType &body_relation, Args &&...args)
			: InteractionDynamics<LocalDynamicsType, DynamicsRange>(body_relation, std::forward<Args>(args)...) {}
		virtual ~InteractionDynamicsWithUpdate(){};

		virtual void exec(Real dt = 0.0) override
		{
			InteractionDynamics<LocalDynamicsType, DynamicsRange>::exec(dt);
			particle_for(
				this->dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->update(i, delta); },
				dt);
		};

		virtual void parallel_exec(Real dt = 0.0) override
		{
			InteractionDynamics<LocalDynamicsType, DynamicsRange>::parallel_exec(dt);
			particle_parallel_for(
				this->dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->update(i, delta); },
				dt);
		};
	};

	/**
	 * @class NewInteractionDynamics1Level
	 * @brief This class includes three steps, including initialization, interaction and update.
	 * It is the most complex particle dynamics type, 
	 * and is typically for computing the main fluid and solid dynamics.
	 */
	template <class LocalDynamicsType, class DynamicsRange = SPHBody>
	class NewInteractionDynamics1Level : public InteractionDynamicsWithUpdate<LocalDynamicsType, DynamicsRange>
	{
	public:
		template <class BodyRelationType, typename... Args>
		NewInteractionDynamics1Level(BodyRelationType &body_relation, Args &&...args)
			: InteractionDynamicsWithUpdate<LocalDynamicsType, DynamicsRange>(body_relation, std::forward<Args>(args)...) {}
		virtual ~NewInteractionDynamics1Level(){};

		virtual void exec(Real dt = 0.0) override
		{
			this->setBodyUpdated();
			this->setupDynamics(dt);

			particle_for(
				this->dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->initialization(i, delta); },
				dt);

			this->runInteractionStep(dt);

			particle_for(
				this->dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->update(i, delta); },
				dt);
		};

		virtual void parallel_exec(Real dt = 0.0) override
		{
			this->setBodyUpdated();
			this->setupDynamics(dt);

			particle_parallel_for(
				this->dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->initialization(i, delta); },
				dt);

			this->parallel_runInteractionStep(dt);

			particle_parallel_for(
				this->dynamics_range_.LoopRange(),
				[&](size_t i, Real delta)
				{ this->update(i, delta); },
				dt);
		};
	};

	/**
	 * @class CombinedLocalDynamics
	 * @brief This is the class for combining multiple local interaction dynamics,
	 * which share the particle loop but are independent from each other,
	 * aiming to increase computing intensity under the data caching environment
	 */
	template <typename... MultipleLocalDynamics>
	class CombinedLocalInteractionDynamics;

	template <>
	class CombinedLocalInteractionDynamics<> : public LocalDynamics
	{
	public:
		template <class BodyRelationType>
		CombinedLocalInteractionDynamics(BodyRelationType &body_relation) : LocalDynamics(body_relation.getDynamicsRange()){};

		void interaction(size_t index_i, Real dt = 0.0){};
	};

	template <class FirstLocalDynamics, class... OtherLocalDynamics>
	class CombinedLocalInteractionDynamics<FirstLocalDynamics, OtherLocalDynamics...> : public LocalDynamics
	{
	protected:
		FirstLocalDynamics first_local_dynamics_;
		CombinedLocalInteractionDynamics<OtherLocalDynamics...> other_local_dynamics_;

	public:
		template <typename BodyRelationType, typename... FirstArgs, typename... OtherArgs>
		CombinedLocalInteractionDynamics(BodyRelationType &body_relation, FirstArgs &&...first_args, OtherArgs &&...other_args)
			: LocalDynamics(body_relation.getDynamicsRange()),
			  first_local_dynamics_(body_relation, std::forward<FirstArgs>(first_args)...),
			  other_local_dynamics_(body_relation, std::forward<OtherArgs>(other_args)...){};

		virtual void setupDynamics(Real dt = 0.0) override
		{
			first_local_dynamics_.setupDynamics(dt);
			other_local_dynamics_.setupDynamics(dt);
		};

		void interaction(size_t index_i, Real dt = 0.0)
		{
			first_local_dynamics_.interaction(index_i, dt);
			other_local_dynamics_.interaction(index_i, dt);
		};
	};

}
#endif // PARTICLE_DYNAMICS_ALGORITHMS_H