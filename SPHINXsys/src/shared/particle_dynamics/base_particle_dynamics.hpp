/**
* @file 	base_particle_dynamics.hpp
* @brief 	This is the implementation of the template class for 3D build
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once

#include "base_particle_dynamics.h"

namespace SPH {
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
		::ParticleDynamicsWithInnerConfigurations(BodyType* body) 
		: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body) {
		current_inner_configuration_ = &body->current_inner_configuration_;
		reference_inner_configuration_ = &body->reference_inner_configuration_;
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsWithContactConfigurations(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>(body), interacting_bodies_(interacting_bodies) {
		/** contact configuration data from the body*/
		SPHBodyVector contact_bodies = body->contact_map_.second;
		ContactParticleList *indexes_contact_particles = &(body->indexes_contact_particles_);
		ContactNeighborList *current_contact_configuration = &(body->current_contact_configuration_);
		ReferenceContactNeighborList *reference_contact_configuration = &(body->reference_contact_configuration_);
		/** finding the interacing bodies */
		for (size_t i = 0; i != interacting_bodies.size(); ++i)
			for (size_t j = 0; j != contact_bodies.size(); ++j) {
				if (static_cast<SPHBody*>(interacting_bodies_[i]) == contact_bodies[j]) {
					interacting_particles_.push_back(dynamic_cast<InteractingParticlesType*>(contact_bodies[j]->base_particles_->PointToThisObject()));
					interacting_material_.push_back(dynamic_cast<InteractingMaterialType*>(contact_bodies[j]->base_material_->PointToThisObject()));
					indexes_interacting_particles_.push_back(&(*indexes_contact_particles)[j]);
					current_interacting_configuration_.push_back(&(*current_contact_configuration)[j]);
					reference_interacting_configuration_.push_back(&(*reference_contact_configuration)[j]);
				}
			}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsByCells<BodyType, ParticlesType, MaterialType>
		::ParticleDynamicsByCells(BodyType* body) 
		: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body)
	{
		mesh_cell_linked_list_ = body->mesh_cell_linked_list_;
		cell_linked_lists_ = mesh_cell_linked_list_->cell_linked_lists_;
		number_of_cells_ = mesh_cell_linked_list_->GetNumberOfCells();
		cell_spacing_ = mesh_cell_linked_list_->GetCellSpacing();
		mesh_lower_bound_ = mesh_cell_linked_list_->GetLowerBound();
		mesh_upper_bound_ = mesh_cell_linked_list_->GetUpperBound();
	}
	//===============================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t number_of_particles, ReturnType temp, 
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &ruduce_operation, Real dt)
	{
		for (size_t i = 0; i < number_of_particles; ++i)
		{
			temp = ruduce_operation(temp, reduce_functor(i, dt));
		}
		return temp;
	}
	//===============================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &ruduce_operation, Real dt)
	{
		temp = parallel_reduce(blocked_range<size_t>(0, number_of_particles),
			temp, [&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType {
			for (size_t i = r.begin(); i != r.end(); ++i) {
				temp0 = ruduce_operation(temp0, reduce_functor(i, dt));
			}
			return temp0;
		},
			[&](ReturnType x, ReturnType y)->ReturnType {
			return ruduce_operation(x, y);
		}
		);
		return temp;
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsSimple<BodyType, ParticlesType, MaterialType>::exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->SetupDynamics(dt);
		InnerIterator(number_of_particles, functor_update_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsSimple<BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->SetupDynamics(dt);
		InnerIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//===============================================================//
	template <class ReturnType, class ReduceOperation, 
		class BodyType, class ParticlesType, class MaterialType>
		ReturnType ParticleDynamicsReduce<ReturnType, ReduceOperation, 
		BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		ReturnType temp = initial_reference_;
		this->SetupReduce();
		temp = ReduceIterator(number_of_particles, 
			temp, functor_reduce_function_, reduce_operation_, dt);
		return this->OutputResult(temp);
	}	
	//===============================================================//
	template <class ReturnType, typename ReduceOperation,
		class BodyType, class ParticlesType, class MaterialType>
		ReturnType ParticleDynamicsReduce<ReturnType, ReduceOperation,
		BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		ReturnType temp = initial_reference_;
		this->SetupReduce();
		temp = ReduceIterator_parallel(number_of_particles,
			temp, functor_reduce_function_, reduce_operation_, dt);
		return this->OutputResult(temp);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>::exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->SetupDynamics(dt);
		InnerIterator(number_of_particles, functor_inner_interaction_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInner<BodyType, ParticlesType, MaterialType>::parallel_exec(Real dt)
	{
		size_t number_of_particles = this->body_->number_of_particles_;
		this->SetupDynamics(dt);
		InnerIterator_parallel(number_of_particles, functor_inner_interaction_, dt);
	}
	//===================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType, MaterialType>::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting(by_cell_lists_particle_indexes_, functor_inner_interaction_, dt);
	}
	//===================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType, MaterialType>::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting_parallel(by_cell_lists_particle_indexes_, functor_inner_interaction_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType, 
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	void ParticleDynamicsContact<BodyType, ParticlesType, MaterialType, 
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		ContactIterator(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsContact<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		ContactIterator_parallel(this->indexes_interacting_particles_, functor_contact_interaction_, dt);
	}
	//===============================================================//
}
