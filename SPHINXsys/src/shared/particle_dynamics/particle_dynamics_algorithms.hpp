/**
* @file 	particle_dynamics_algorithms.hpp
* @brief 	This is the implementation of the template class particle dynamics algorithms
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once
#include "particle_dynamics_algorithms.h"
//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsSimple<BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->setupDynamics(dt);
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsSimple<BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->setupDynamics(dt);
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class ReturnType, class ReduceOperation, class BodyType, class ParticlesType, class MaterialType>
	ReturnType ParticleDynamicsReduce<ReturnType, ReduceOperation, BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->SetupReduce();
		ReturnType temp = ReduceIterator(number_of_particles,
			initial_reference_, functor_reduce_function_, reduce_operation_, dt);
		return this->OutputResult(temp);
	}
	//=================================================================================================//
	template <class ReturnType, typename ReduceOperation, class BodyType, class ParticlesType, class MaterialType>
	ReturnType ParticleDynamicsReduce<ReturnType, ReduceOperation, BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->SetupReduce();
		ReturnType temp = ReduceIterator_parallel(number_of_particles,
			initial_reference_, functor_reduce_function_, reduce_operation_, dt);
		return this->OutputResult(temp);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->setupDynamics(dt);
		InnerIterator(number_of_particles, functor_inner_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->setupDynamics(dt);
		InnerIterator_parallel(number_of_particles, functor_inner_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsInnerWithUpdate<BodyType, ParticlesType, MaterialType>
	::ParticleDynamicsInnerWithUpdate(BodyType* body)
	: ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>(body),
		functor_update_(std::bind(&ParticleDynamicsInnerWithUpdate::Update, this, _1, _2)) {}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerWithUpdate<BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>::exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerWithUpdate<BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>::parallel_exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsInner1Level<BodyType, ParticlesType, MaterialType>
	::ParticleDynamicsInner1Level(BodyType* body)
	: ParticleDynamicsInnerWithUpdate<BodyType, ParticlesType, MaterialType>(body),
		functor_initialization_(std::bind(&ParticleDynamicsInner1Level::Initialization, this, _1, _2)) {}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner1Level<BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_initialization_, dt);
		InnerIterator(number_of_particles, this->functor_inner_interaction_, dt);
		InnerIterator(number_of_particles, this->functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner1Level<BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_inner_interaction_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsContact<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsContact<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	::ParticleDynamicsComplex(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
	: ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
		functor_complex_interaction_(std::bind(&ParticleDynamicsComplex::ComplexInteraction, this, _1, _2)) 
	{}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_complex_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_complex_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	::ParticleDynamicsComplexWithUpdate(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies)
	: ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
		functor_update_(std::bind(&ParticleDynamicsComplexWithUpdate::Update, this, _1, _2)) 
	{}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>::exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
			InteractingBodyType, InteractingParticlesType, InteractingMaterialType>::parallel_exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
	::ParticleDynamicsComplex1Level(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies)
	: ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies),
		functor_initialization_(std::bind(&ParticleDynamicsComplex1Level::Initialization, this, _1, _2)) 
	{}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_initialization_, dt);
		InnerIterator(number_of_particles, this->functor_complex_interaction_, dt);
		InnerIterator(number_of_particles, this->functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_complex_interaction_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_update_, dt);
	}
	//===================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		ParticleDynamicsComplexSplit<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsComplexSplit(BodyType* body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>(body, interacting_bodies) {}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexSplit<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, this->functor_initialization_, dt);
		InnerIteratorSplitting(this->split_cell_lists_, this->functor_complex_interaction_, dt);
		InnerIterator(number_of_particles, this->functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexSplit<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, this->functor_initialization_, dt);
		InnerIteratorSplitting_parallel(this->split_cell_lists_, this->functor_complex_interaction_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_update_, dt);
	}
	//===============================================================//
}
//=================================================================================================//