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

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH
{
	/**
	* @class ParticleDynamicsSimple
	* @brief Simple particle dynamics without considering particle interaction
	*/
	class ParticleDynamicsSimple : public ParticleDynamics<void>
	{
	public:
		explicit ParticleDynamicsSimple(SPHBody &sph_body)
			: ParticleDynamics<void>(sph_body),
			  functor_update_(std::bind(&ParticleDynamicsSimple::Update, this, _1, _2)){};
		virtual ~ParticleDynamicsSimple(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_update_;
	};

	/**
	* @class ParticleDynamicsReduce
	* @brief Base abstract class for reduce
	*/
	template <class ReturnType, typename ReduceOperation>
	class ParticleDynamicsReduce : public ParticleDynamics<ReturnType>
	{
	public:
		explicit ParticleDynamicsReduce(SPHBody &sph_body)
			: ParticleDynamics<ReturnType>(sph_body), quantity_name_("ReducedQuantity"), initial_reference_(),
			  functor_reduce_function_(std::bind(&ParticleDynamicsReduce::ReduceFunction, this, _1, _2)){};
		virtual ~ParticleDynamicsReduce(){};

		ReturnType InitialReference() { return initial_reference_; };
		std::string QuantityName() { return quantity_name_; };

		virtual ReturnType exec(Real dt = 0.0) override
		{
			size_t total_real_particles = this->base_particles_->total_real_particles_;
			this->setBodyUpdated();
			SetupReduce();
			ReturnType temp = ReduceIterator(total_real_particles,
											 initial_reference_, functor_reduce_function_, reduce_operation_, dt);
			return OutputResult(temp);
		};
		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			size_t total_real_particles = this->base_particles_->total_real_particles_;
			this->setBodyUpdated();
			SetupReduce();
			ReturnType temp = ReduceIterator_parallel(total_real_particles,
													  initial_reference_, functor_reduce_function_, reduce_operation_, dt);
			return this->OutputResult(temp);
		};

	protected:
		ReduceOperation reduce_operation_;
		std::string quantity_name_;

		/** initial or reference value */
		ReturnType initial_reference_;
		virtual void SetupReduce(){};
		virtual ReturnType ReduceFunction(size_t index_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
		ReduceFunctor<ReturnType> functor_reduce_function_;
	};

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
		StdVec<ParticleDynamics<void> *> pre_processes_;
		/** post process such as impose constraint */
		StdVec<ParticleDynamics<void> *> post_processes_;

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
		SplitCellLists &split_cell_lists_;
	};

	//----------------------------------------------------------------------
	//		New version particle dynamics base classes.
	//		Aiming to replace old version in the future.
	//		In this version, std functors and dynamic bindings are replaced by
	//		standard functors with static binding
	//----------------------------------------------------------------------

	/** loop for particle dynamics. sequential computing. */
	template <typename LocalDynamicsType>
	void particle_for(size_t total_real_particles, LocalDynamicsType &local_dynamics, Real dt = 0.0)
	{
		for (size_t i = 0; i < total_real_particles; ++i)
			local_dynamics(i, dt);
	};
	/** loop for particle dynamics.  parallel computing. */
	template <typename LocalDynamicsType>
	void particle_parallel_for(size_t total_real_particles, LocalDynamicsType &local_dynamics, Real dt = 0.0)
	{
		parallel_for(
			blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					local_dynamics(i, dt);
				}
			},
			ap);
	};

	template <class LocalDynamicsType>
	class SimpleDynamics : public ParticleDynamics<void>
	{
	public:
		LocalDynamicsType local_dynamics_;

		explicit SimpleDynamics(SPHBody &sph_body)
			: ParticleDynamics<void>(sph_body), local_dynamics_(sph_body){};
		virtual ~SimpleDynamics(){};

		virtual void exec(Real dt = 0.0) override
		{
			setBodyUpdated();
			setupDynamics(dt);
			size_t total_real_particles = base_particles_->total_real_particles_;
			particle_for(total_real_particles, local_dynamics_, dt);
		};

		virtual void parallel_exec(Real dt = 0.0) override
		{
			setBodyUpdated();
			setupDynamics(dt);
			size_t total_real_particles = base_particles_->total_real_particles_;
			particle_parallel_for(total_real_particles, local_dynamics_, dt);
		};
	};

	/** loop particles for reducing. sequential computing. */
	template <class ReturnType, typename ReduceOperation, class LocalEvaluation>
	ReturnType particle_reduce(size_t total_real_particles, ReturnType temp,
							   ReduceOperation &reduce_operation, LocalEvaluation &local_evaluation, Real dt = 0.0)
	{
		for (size_t i = 0; i < total_real_particles; ++i)
		{
			temp = reduce_operation(temp, local_evaluation(i, dt));
		}
		return temp;
	};

	/** loop particles for reducing. parallel computing. */
	template <class ReturnType, typename ReduceOperation, class LocalEvaluation>
	ReturnType particle_parallel_reduce(size_t total_real_particles, ReturnType temp,
										ReduceOperation &reduce_operation, LocalEvaluation &local_evaluation, Real dt = 0.0)
	{
		return parallel_reduce(
			blocked_range<size_t>(0, total_real_particles),
			temp, [&](const blocked_range<size_t> &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = reduce_operation(temp0, local_evaluation(i, dt));
				}
				return temp0;
			},
			[&](ReturnType x, ReturnType y) -> ReturnType
			{
				return reduce_operation(x, y);
			});
	};

	/**
	* @class BaseReduce
	* @brief Base class for reduce operations
	* @details Specific reduce methods will be derived from this class.
	*		However, since I try to not use virtual functions for local particle based operations,
	*		there is no virtual function to be overriden.
	*/
	template <typename ValueType, typename ReduceOperationType>
	class BaseReduce
	{
	public:
		std::string quantity_name_;
		ValueType initial_reference_;
		ReduceOperationType reduce_operation_;
		typedef ValueType ReducedValueType;

		BaseReduce(const std::string &quantity_name, const ValueType &initial_reference)
			: quantity_name_(quantity_name), initial_reference_(initial_reference){};
		virtual ~BaseReduce(){};
	};

	/**
	* @class ReduceDynamics
	* @brief Template class for reduce operation through all particles
	* @details For this, the template LocalReducedMethod should be defined with signature functions
	*		and operater overloading.
	*/
	template <class LocalReduceMethod, typename ReturnType = typename LocalReduceMethod::ReducedValueType>
	class ReduceDynamics : public ParticleDynamics<ReturnType>
	{
	public:
		LocalReduceMethod local_reduce_;

		explicit ReduceDynamics(SPHBody &sph_body)
			: ParticleDynamics<ReturnType>(sph_body), local_reduce_(sph_body),
			  initial_reference_(local_reduce_.initial_reference_),
			  quantity_name_(local_reduce_.quantity_name_){};
		virtual ~ReduceDynamics(){};

		ReturnType InitialReference() { return initial_reference_; };
		std::string QuantityName() { return quantity_name_; };

		virtual ReturnType exec(Real dt = 0.0) override
		{
			size_t total_real_particles = this->base_particles_->total_real_particles_;
			local_reduce_.setupReduce();
			ReturnType result = particle_reduce(total_real_particles, initial_reference_,
												local_reduce_, local_reduce_.operation_, dt);
			return local_reduce_.outputResult(result);
		};
		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			size_t total_real_particles = this->base_particles_->total_real_particles_;
			local_reduce_.setupReduce();
			ReturnType result = particle_parallel_reduce(total_real_particles, initial_reference_,
														 local_reduce_, local_reduce_.operation_, dt);
			return local_reduce_.outputResult(result);
		};

	protected:
		/** inital or reference value */
		ReturnType initial_reference_;
		std::string quantity_name_;
	};
}
#endif //PARTICLE_DYNAMICS_ALGORITHMS_H