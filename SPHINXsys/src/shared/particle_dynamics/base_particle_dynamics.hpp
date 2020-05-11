/**
* @file 	base_particle_dynamics.hpp
* @brief 	This is the implementation of the template class for 3D build
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once

#include "base_particle_dynamics.h"
//=================================================================================================//
namespace SPH {
//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>
		::ParticleDynamicsWithInnerConfigurations(BodyType* body) 
		: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body) {
		inner_configuration_ = &body->inner_configuration_;
	}
//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
	ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::ParticleDynamicsWithContactConfigurations(BodyType *body, StdVec<InteractingBodyType*> interacting_bodies)
		: ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType, MaterialType>(body), interacting_bodies_(interacting_bodies) {
		/** contact configuration data from the body*/
		SPHBodyVector contact_bodies = body->contact_map_.second;
		ContactParticles& indexes_contact_particles = body->indexes_contact_particles_;
		ContatcParticleConfiguration& current_contact_configuration = body->contact_configuration_;
		/** finding the interacing bodies */
		for (size_t i = 0; i != interacting_bodies.size(); ++i)
			for (size_t j = 0; j != contact_bodies.size(); ++j) {
				if (static_cast<SPHBody*>(interacting_bodies_[i]) == contact_bodies[j]) {
					interacting_particles_.push_back(dynamic_cast<InteractingParticlesType*>(contact_bodies[j]->base_particles_->PointToThisObject()));
					interacting_material_.push_back(dynamic_cast<InteractingMaterialType*>(contact_bodies[j]->base_particles_->base_material_->PointToThisObject()));
					indexes_interacting_particles_.push_back(&(indexes_contact_particles[j]));
					current_interacting_configuration_.push_back(&(current_contact_configuration[j]));
				}
			}
	}
//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsByCells<BodyType, ParticlesType, MaterialType>
		::ParticleDynamicsByCells(BodyType* body) 
		: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body)
	{
		mesh_cell_linked_list_ = body->base_mesh_cell_linked_list_;
		cell_linked_lists_ = mesh_cell_linked_list_->getCellLinkedLists();
		number_of_cells_ = mesh_cell_linked_list_->getNumberOfCells();
		cell_spacing_ = mesh_cell_linked_list_->getCellSpacing();
		mesh_lower_bound_ = mesh_cell_linked_list_->getMeshLowerBound();
	}
//=================================================================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &reduce_operation, Real dt)
	{
		for (size_t i = 0; i < number_of_particles; ++i)
		{
			temp = reduce_operation(temp, reduce_functor(i, dt));
		}
		return temp;
	}
//=================================================================================================//
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &reduce_operation, Real dt)
	{
		return parallel_reduce(blocked_range<size_t>(0, number_of_particles),
			temp, [&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType {
			for (size_t i = r.begin(); i != r.end(); ++i) {
				temp0 = reduce_operation(temp0, reduce_functor(i, dt));
			}
			return temp0;
		},
			[&](ReturnType x, ReturnType y)->ReturnType {
			return reduce_operation(x, y);
		}
		);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	ParticleDynamicsCellListSplitting<BodyType, ParticlesType, MaterialType>
		::ParticleDynamicsCellListSplitting(BodyType* body)
		: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body),
		functor_cell_list_(std::bind(&ParticleDynamicsCellListSplitting::CellListInteraction, this, _1, _2))
	{
		mesh_cell_linked_list_ = body->base_mesh_cell_linked_list_;
		cell_linked_lists_ = mesh_cell_linked_list_->getCellLinkedLists();
		number_of_cells_ = mesh_cell_linked_list_->getNumberOfCells();
		cell_spacing_ = mesh_cell_linked_list_->getCellSpacing();
		mesh_lower_bound_ = mesh_cell_linked_list_->getMeshLowerBound();
		cutoff_radius_ = mesh_cell_linked_list_->getCellSpacing();
		kernel_ = body->kernel_;
	};
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsCellListSplitting<BodyType, ParticlesType, MaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		CellListIteratorSplitting(this->split_cell_lists_, functor_cell_list_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsCellListSplitting<BodyType, ParticlesType, MaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		CellListIteratorSplitting_parallel(this->split_cell_lists_, functor_cell_list_, dt);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType, MaterialType>::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting(this->split_cell_lists_, functor_inner_interaction_, dt);
	}
//=================================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType, MaterialType>::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting_parallel(this->split_cell_lists_, functor_inner_interaction_, dt);
	}
	//=============================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting(this->split_cell_lists_, functor_particle_interaction_, dt);
	}
	//=============================================================================================//
	template <class BodyType, class ParticlesType, class MaterialType,
		class InteractingBodyType, class InteractingParticlesType, class InteractingMaterialType>
		void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, MaterialType,
		InteractingBodyType, InteractingParticlesType, InteractingMaterialType>
		::parallel_exec(Real dt)
	{
		this->SetupDynamics(dt);
		InnerIteratorSplitting_parallel(this->split_cell_lists_, functor_particle_interaction_, dt);
	}
}
//=================================================================================================//
