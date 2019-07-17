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
	* @class ParticleDynamicsInner1Level
	* @brief This class includes an initialization, an inner interaction and an update steps
	*/
	template <class BodyType, class ParticlesType, class MaterialType = Material>
	class ParticleDynamicsInner1Level 
		: public ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
	{
	protected:
		size_t number_of_particles_;

		virtual void Initialization(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_initialization_;
		InnerFunctor functor_inner_interaction_;
		InnerFunctor functor_update_;
	public:
		ParticleDynamicsInner1Level(BodyType* body);
		virtual ~ParticleDynamicsInner1Level() {};

		virtual void exec(Real dt = 0.0);
		virtual void parallel_exec(Real dt = 0.0);
	};

	/**
	* @class ParticleDynamicsComplex
	* @brief compex operations involving both inner and contact particle dynamics
	*/
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = Material>
	class ParticleDynamicsComplex 
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		size_t number_of_particles_;

		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		InnerFunctor functor_inner_interaction_;
		ContactFunctor functor_contact_interaction_;
	public:
		ParticleDynamicsComplex(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplex() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplexWithUpdate
	* @brief an extra particle updating is inlcuded
	*/
	template <class BodyType, class ParticlesType, class MaterialType, 
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = Material>
	class ParticleDynamicsComplexWithUpdate 
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		size_t number_of_particles_;

		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_inner_interaction_;
		ContactFunctor functor_contact_interaction_;
		InnerFunctor functor_update_;
	public:
		ParticleDynamicsComplexWithUpdate(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplexWithUpdate() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplex1Level
	* @brief compex operations involving one level with both inner
	* and contact particle dynamics, initilaize and update
	*/
	template <class BodyType, class ParticlesType, class MaterialType, 
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = Material>
	class ParticleDynamicsComplex1Level 
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		size_t number_of_particles_;

		virtual void Initialization(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_initialization_;
		InnerFunctor functor_inner_interaction_;
		ContactFunctor functor_contact_interaction_;
		InnerFunctor functor_update_;

	public:
		ParticleDynamicsComplex1Level(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplex1Level() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplex2Levels
	* @brief compex operations involving 2 level with both inner
	* and contact particle dynamics, initilaize and update
	*/
	template <class BodyType, class ParticlesType, class MaterialType, 
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = Material>
	class ParticleDynamicsComplex2Levels 
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		size_t number_of_particles_;

		virtual void Initialization(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		virtual void Intermediate(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void InnerInteraction2nd(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction2nd(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		InnerFunctor functor_initialization_;
		InnerFunctor functor_inner_interaction_;
		ContactFunctor functor_contact_interaction_;
		InnerFunctor functor_intermediate_;
		InnerFunctor functor_inner_interaction_2nd_;
		ContactFunctor functor_contact_interaction_2nd_;
		InnerFunctor functor_update_;
	public:
		ParticleDynamicsComplex2Levels(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplex2Levels() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class ParticleDynamicsComplexSplitting
	* @brief inner interaction use splitting scheme
	*/
	template <class BodyType, class ParticlesType, class MaterialType, 
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType = Material>
	class ParticleDynamicsComplexSplitting 
		: public ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	{
	protected:
		ByCellLists by_cell_lists_particle_indexes_;

		virtual void SetupDynamics(Real dt = 0.0) {};
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual void ContactInteraction(size_t index_particle_i, size_t interacting_body_index, Real dt = 0.0) = 0;
		InnerFunctor functor_inner_interaction_;
		ContactFunctor functor_contact_interaction_;
	public:
		explicit ParticleDynamicsComplexSplitting(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies);
		virtual ~ParticleDynamicsComplexSplitting() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
}

