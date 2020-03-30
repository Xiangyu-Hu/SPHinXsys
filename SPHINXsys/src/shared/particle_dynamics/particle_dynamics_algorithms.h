/**
* @file 	particle_dynamics_algorithms.h
* @brief 	This is the classes for algorithms which are using at least two particle dynamics
* and some generic algorithms used for particle operations
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH {

	/**
	* @class ParticleDynamicsSimple
	* @brief Simple particle dynamics base class
	*/
	template <class BodyType, class ParticlesType = BaseParticles, class MaterialType = BaseMaterial>
	class ParticleDynamicsSimple : public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_update_;
	public:
		explicit ParticleDynamicsSimple(BodyType* body)
			: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body),
			functor_update_(std::bind(&ParticleDynamicsSimple::Update, this, _1, _2)) {};
		virtual ~ParticleDynamicsSimple() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsReduce
	* @brief Base abstract class for reduce
	*/
	template <class ReturnType, typename ReduceOperation,
		class BodyType, class ParticlesType = BaseParticles, class MaterialType = BaseMaterial>
		class ParticleDynamicsReduce : public ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>
	{
	protected:
		ReduceOperation reduce_operation_;

		/** inital or refence value */
		ReturnType initial_reference_;
		virtual void SetupReduce() {};
		virtual ReturnType ReduceFunction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
		ReduceFunctor<ReturnType> functor_reduce_function_;
	public:
		explicit ParticleDynamicsReduce(BodyType* body)
			: ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>(body),
			functor_reduce_function_(std::bind(&ParticleDynamicsReduce::ReduceFunction, this, _1, _2)),
			initial_reference_() {};
		virtual ~ParticleDynamicsReduce() {};

		virtual ReturnType exec(Real dt = 0.0) override;
		virtual ReturnType parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsInner
	* @brief This is the class for inner interactions
	* in which one the particles from the same body
	* interact with each other
	*/
	template <class BodyType, class ParticlesType = BaseParticles, class MaterialType = BaseMaterial>
	class ParticleDynamicsInner : public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
	{
	protected:
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_inner_interaction_;
	public:
		explicit ParticleDynamicsInner(BodyType* body) :
			ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>(body),
			functor_inner_interaction_(std::bind(&ParticleDynamicsInner::InnerInteraction, this, _1, _2)) {};
		virtual ~ParticleDynamicsInner() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsInnerWithUpdate
	* @brief This class includes an initialization, an inner interaction and an update steps
	*/
	template <class BodyType, class ParticlesType, class MaterialType = BaseMaterial>
	class ParticleDynamicsInnerWithUpdate
		: public ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_update_;
	public:
		ParticleDynamicsInnerWithUpdate(BodyType* body);
		virtual ~ParticleDynamicsInnerWithUpdate() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsInner1Level
	* @brief This class includes an initialization, an inner interaction and an update steps
	*/
	template <class BodyType, class ParticlesType, class MaterialType = BaseMaterial>
	class ParticleDynamicsInner1Level 
		: public ParticleDynamicsInnerWithUpdate<BodyType, ParticlesType, MaterialType>
	{
	protected:
		virtual void Initialization(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_initialization_;
	public:
		ParticleDynamicsInner1Level(BodyType* body);
		virtual ~ParticleDynamicsInner1Level() {};

		virtual void exec(Real dt = 0.0);
		virtual void parallel_exec(Real dt = 0.0);
	};

	/**
	 * @class ParticleDynamicsContact
	 * @brief This is the class for contact interactions
	 */
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = BaseMaterial>
		class ParticleDynamicsContact
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		ContactFunctor functor_contact_interaction_;
	public:
		explicit ParticleDynamicsContact(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies)
			: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
			functor_contact_interaction_(std::bind(&ParticleDynamicsContact::ContactInteraction, this, _1, _2, _3)) {};
		virtual ~ParticleDynamicsContact() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplex
	* @brief compex operations combining both inner and contact particle dynamics together
	*/
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = BaseMaterial>
		class ParticleDynamicsComplex
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_complex_interaction_;

	public:
		ParticleDynamicsComplex(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplex() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplexWithUpdate
	* @brief compex operations combining both inner and contact particle dynamics together
	*/
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = BaseMaterial>
		class ParticleDynamicsComplexWithUpdate
		: public ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_update_;

	public:
		ParticleDynamicsComplexWithUpdate(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplexWithUpdate() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplex1Level
	* @brief compex operations combining both inner and contact particle dynamics together
	*/
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = BaseMaterial>
		class ParticleDynamicsComplex1Level
		: public ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		virtual void Initialization(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_initialization_;

	public:
		ParticleDynamicsComplex1Level(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplex1Level() {};
	
		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplexSplit
	* @brief compex split operations combining both inner and contact particle dynamics together
	*/
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = BaseMaterial>
		class ParticleDynamicsComplexSplit
		: public ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	public:
		ParticleDynamicsComplexSplit(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplexSplit() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
}
