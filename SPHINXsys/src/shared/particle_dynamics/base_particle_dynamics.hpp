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
	template <class BodyType, class ParticlesType>
	ParticleDynamicsWithInnerConfigurations<BodyType, ParticlesType>
		::ParticleDynamicsWithInnerConfigurations(BodyType* body) : ParticleDynamics(body) {
		current_inner_configuration_ = &body->current_inner_configuration_;
		reference_inner_configuration_ = &body->reference_inner_configuration_;
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	ParticleDynamicsWithContactConfigurations<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::ParticleDynamicsWithContactConfigurations(BodyType *body, StdVec<InteractingBodytype*> interacting_bodies)
		: ParticleDynamicsWithInnerConfigurations(body), interacting_bodies_(interacting_bodies) {
		/** contact configuration data from the body*/
		SPHBodyVector contact_bodies = body->contact_map_.second;
		ContactParticleList *indexes_contact_particles = &(body->indexes_contact_particles_);
		ContactNeighborList *current_contact_configuration = &(body->current_contact_configuration_);
		ReferenceContactNeighborList *reference_contact_configuration = &(body->reference_contact_configuration_);
		/** finding the interacing bodies */
		for (size_t i = 0; i != interacting_bodies.size(); ++i)
			for (size_t j = 0; j != contact_bodies.size(); ++j) {
				if (static_cast<SPHBody*>(interacting_bodies_[i]) == contact_bodies[j]) {
					interacting_particles_.push_back(dynamic_cast<InteractingParticlesType*>(contact_bodies[j]->base_particles_.PointToThisObject()));
					indexes_interacting_particles_.push_back(&(*indexes_contact_particles)[j]);
					current_interacting_configuration_.push_back(&(*current_contact_configuration)[j]);
					reference_interacting_configuration_.push_back(&(*reference_contact_configuration)[j]);
				}
			}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType>
	ParticleDynamicsByCells<BodyType, ParticlesType>
		::ParticleDynamicsByCells(BodyType* body) : ParticleDynamics<void, BodyType, ParticlesType>(body),
		mesh_cell_linked_list_(*(body->mesh_cell_linked_list_))
	{
		cell_linked_lists_ = mesh_cell_linked_list_.cell_linked_lists_;
		number_of_cells_ = mesh_cell_linked_list_.GetNumberOfCells();
		cell_spacing_ = mesh_cell_linked_list_.GetCellSpacing();
		mesh_lower_bound_ = mesh_cell_linked_list_.GetLowerBound();
		mesh_upper_bound_ = mesh_cell_linked_list_.GetUpperBound();
	}
	//===============================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsSimple<BodyType, ParticlesType>::exec(Real dt)
	{
		SetupDynamics(dt);
		for (size_t i = 0; i < number_of_particles_; ++i) ParticleUpdate(i, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsSimple<BodyType, ParticlesType>::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) ParticleUpdate(i, dt);
		}, ap);
	}
	//===============================================================//
	template <class ReturnType, class BodyType, class ParticlesType>
	ReturnType ParticleDynamicsReduce<ReturnType, BodyType, ParticlesType>::exec(Real dt) 
	{
		ReturnType temp = initial_reference_;
		SetupReduce();
		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			temp = ReduceOperation(temp, ReduceFunction(i, dt));
		}
		return OutputResult(temp);
	}	
	//===============================================================//
	template <class ReturnType, class BodyType, class ParticlesType>
	ReturnType ParticleDynamicsReduce<ReturnType, BodyType, ParticlesType>::parallel_exec(Real dt) 
	{
		ReturnType temp = initial_reference_;
		SetupReduce();
		temp = parallel_reduce(blocked_range<size_t>(0, number_of_particles_),
			temp, [&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType {
			for (size_t i = r.begin(); i != r.end(); ++i) {
				temp0 = ReduceOperation(temp0, ReduceFunction(i, dt));
			}
			return temp0;
		},
			[this](ReturnType x, ReturnType y)->ReturnType {
			return ReduceOperation(x, y);
		}
		);
		return OutputResult(temp);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsInner<BodyType, ParticlesType>::exec(Real dt)
	{
		SetupDynamics(dt);
		for (size_t i = 0; i < number_of_particles_; ++i) InnerInteraction(i, dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsInner<BodyType, ParticlesType>::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) InnerInteraction(i, dt);
		}, ap);
	}
	//===================================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType>::exec(Real dt)
	{
		SetupDynamics(dt);
		//forward sweeping
		for (size_t k = 0; k < numer_of_lists_; ++k) {
			StdVec<IndexVector> &lists_particle_indexes_
				= body_->by_cell_lists_particle_indexes_[k];
			for (size_t l = 0; l < lists_particle_indexes_.size(); ++l)
			{
				for (size_t i = 0; i < lists_particle_indexes_[l].size(); ++i)
				{
					InnerInteraction(lists_particle_indexes_[l][i], dt);
				}
			}
		}

		//backward sweeping
		for (size_t k = numer_of_lists_ - 1; k >= 0; --k) {
			StdVec<IndexVector> &lists_particle_indexes_
				= body_->by_cell_lists_particle_indexes_[k];
			for (size_t l = lists_particle_indexes_.size() - 1; l >= 0; --l)
			{
				for (size_t i = lists_particle_indexes_[l].size() - 1; i >= 0; --i)
				{
					InnerInteraction(lists_particle_indexes_[l][i], dt);
				}
			}
		}
	}
	//===================================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsInnerSplitting<BodyType, ParticlesType>::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		//forward sweeping
		for (size_t k = 0; k < numer_of_lists_; ++k) {
			StdVec<IndexVector> &lists_particle_indexes_
				= body_->by_cell_lists_particle_indexes_[k];
			parallel_for(blocked_range<size_t>(0, lists_particle_indexes_.size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.begin(); l < r.end(); ++l) {
					for (size_t i = 0; i < lists_particle_indexes_[l].size(); ++i)
					{
						InnerInteraction(lists_particle_indexes_[l][i], dt);
					}
				}
			}, ap);
		}

		//backward sweeping
		for (size_t k = numer_of_lists_; k >= 1; --k) {
			StdVec<IndexVector> &lists_particle_indexes_
				= body_->by_cell_lists_particle_indexes_[k - 1];
			parallel_for(blocked_range<size_t>(0, lists_particle_indexes_.size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.end(); l >= r.begin() + 1; --l) {
					IndexVector &particle_indexes = lists_particle_indexes_[l - 1];
					for (size_t i = particle_indexes.size(); i >= 1; --i)
					{
						InnerInteraction(particle_indexes[i - 1], dt);
					}
				}
			}, ap);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsContact<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>::exec(Real dt)
	{
		SetupDynamics(dt);
		for (size_t k = 0; k < interacting_bodies_.size(); ++k)
			for (size_t l = 0; l < (*indexes_interacting_particles_[k]).size(); ++l) {
				size_t particle_index_i = (*indexes_interacting_particles_[k])[l];
				ContactInteraction(particle_index_i, k, dt);
			}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsContact<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		for (size_t k = 0; k < interacting_bodies_.size(); ++k)
			parallel_for(blocked_range<size_t>(0, (*indexes_interacting_particles_[k]).size()),
				[&](const blocked_range<size_t>& r) {
			for (size_t l = r.begin(); l < r.end(); ++l) {
				size_t particle_index_i = (*indexes_interacting_particles_[k])[l];
				ContactInteraction(particle_index_i, k, dt);
			}
		}, ap);
	}
	//===============================================================//
}
