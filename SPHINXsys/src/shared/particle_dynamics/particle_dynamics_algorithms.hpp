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
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsContact<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_contact_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsContact<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_contact_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_complex_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_complex_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::exec(Real dt)
	{
		ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
			ContactBodyType, ContactParticlesType, ContactMaterialType>::exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::parallel_exec(Real dt)
	{
		ParticleDynamicsComplex<BodyType, ParticlesType, MaterialType,
			ContactBodyType, ContactParticlesType, ContactMaterialType>::parallel_exec(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
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
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, functor_initialization_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_complex_interaction_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
		void ParticleDynamicsComplexSplit<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
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
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
		void ParticleDynamicsComplexSplit<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		size_t number_of_particles = this->body_->number_of_particles_;
		InnerIterator_parallel(number_of_particles, this->functor_initialization_, dt);
		InnerIteratorSplitting_parallel(this->split_cell_lists_, this->functor_complex_interaction_, dt);
		InnerIterator_parallel(number_of_particles, this->functor_update_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsCellListSplitting<BodyType, ParticlesType, MaterialType>
		::ParticleDynamicsCellListSplitting(SPHBody* body)
		: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body),
		functor_cell_list_(std::bind(&ParticleDynamicsCellListSplitting::CellListInteraction, this, _1, _2))
	{
		cutoff_radius_ = this->mesh_cell_linked_list_->CellSpacing();
		kernel_ = body->kernel_;
	};
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsCellListSplitting<BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		CellListIteratorSplitting(this->split_cell_lists_, functor_cell_list_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsCellListSplitting<BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		CellListIteratorSplitting_parallel(this->split_cell_lists_, functor_cell_list_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType, MaterialType>::exec(Real dt)
	{
		this->setupDynamics(dt);
		InnerIteratorSplittingSweeping(this->split_cell_lists_, functor_inner_interaction_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType, MaterialType>::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		InnerIteratorSplittingSweeping_parallel(this->split_cell_lists_, functor_inner_interaction_, dt);
	}
	//=============================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
		void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		InnerIteratorSplittingSweeping(this->split_cell_lists_, functor_particle_interaction_, dt);
	}
	//=============================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class ContactBodyType, class ContactParticlesType, class ContactMaterialType>
		void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, MaterialType,
		ContactBodyType, ContactParticlesType, ContactMaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		InnerIteratorSplittingSweeping_parallel(this->split_cell_lists_, functor_particle_interaction_, dt);
	}
	//=============================================================================================//
}
//=================================================================================================//