/**
* @file 	particle_dynamics_algorithms.hpp
* @brief 	This is the implementation of the template class particle dynamics algorithms
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "particle_dynamics_algorithms.h"

namespace SPH{
	//===================================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsInner1Level<BodyType, ParticlesType>::exec(Real dt)
	{
		SetupDynamics(dt);
		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			Initialization(i, dt);
		}

		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			InnerInteraction(i, dt);
		}

		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			Update(i, dt);
		}
	}
	//===================================================================//
	template <class BodyType, class ParticlesType>
	void ParticleDynamicsInner1Level<BodyType, ParticlesType>::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Initialization(i, dt);
			}
		}, ap);

		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				InnerInteraction(i, dt);
			}
		}, ap);

		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Update(i, dt);
			}
		}, ap);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInnerInteraction_exec(Real dt)
	{
		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			InnerInteraction(i, dt);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingContactInteraction_exec(Real dt)
	{
		for (size_t k = 0; k < interacting_bodies_.size(); ++k)
		{
			for (size_t l = 0; l < (*indexes_interacting_particles_[k]).size(); ++l)
			{
				size_t particle_index_i = (*indexes_interacting_particles_[k])[l];
				ContactInteraction(particle_index_i, k, dt);
			}
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::exec(Real dt)
	{
		SetupDynamics(dt);
		LoopingInnerInteraction_exec(dt);
		LoopingContactInteraction_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInnerInteraction_parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				InnerInteraction(i, dt);
			}
		}, ap);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingContactInteraction_parallel_exec(Real dt)
	{
		for (size_t k = 0; k < interacting_bodies_.size(); ++k)
		{

			parallel_for(blocked_range<size_t>(0, (*indexes_interacting_particles_[k]).size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.begin(); l < r.end(); ++l) {
					size_t particle_index_i = (*indexes_interacting_particles_[k])[l];
					ContactInteraction(particle_index_i, k, dt);
				}
			}, ap);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		LoopingInnerInteraction_parallel_exec(dt);
		LoopingContactInteraction_parallel_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingUpdate_exec(Real dt)
	{
		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			Update(i, dt);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingUpdate_parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Update(i, dt);
			}
		}, ap);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::exec(Real dt)
	{
		SetupDynamics(dt);
		LoopingInnerInteraction_exec(dt);
		LoopingContactInteraction_exec(dt);
		LoopingUpdate_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplexWithUpdate<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		LoopingInnerInteraction_parallel_exec(dt);
		LoopingContactInteraction_parallel_exec(dt);
		LoopingUpdate_parallel_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInitialization_exec(Real dt)
	{
		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			Initialization(i, dt);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::exec(Real dt)
	{

		SetupDynamics(dt);
		LoopingInitialization_exec(dt);
		LoopingInnerInteraction_exec(dt);
		LoopingContactInteraction_exec(dt);
		LoopingUpdate_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInitialization_parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Initialization(i, dt);
			}
		}, ap);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex1Level<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		LoopingInitialization_parallel_exec(dt);
		LoopingInnerInteraction_parallel_exec(dt);
		LoopingContactInteraction_parallel_exec(dt);
		LoopingUpdate_parallel_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingIntermediate_exec(Real dt)
	{
		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			Intermediate(i, dt);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInnerInteraction2nd_exec(Real dt)
	{
		for (size_t i = 0; i < number_of_particles_; ++i)
		{
			InnerInteraction2nd(i, dt);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingContactInteraction2nd_exec(Real dt)
	{
		for (size_t k = 0; k < interacting_bodies_.size(); ++k)
		{
			for (size_t l = 0; l < (*indexes_interacting_particles_[k]).size(); ++l)
			{
				size_t particle_index_i = (*indexes_interacting_particles_[k])[l];
				ContactInteraction2nd(particle_index_i, k, dt);
			}
		}
	}	
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::exec(Real dt)
	{
		SetupDynamics(dt);
		LoopingInitialization_exec(dt);
		LoopingInnerInteraction_exec(dt);
		LoopingContactInteraction_exec(dt);
		LoopingIntermediate_exec(dt);
		LoopingInnerInteraction2nd_exec(dt);
		LoopingContactInteraction2nd_exec(dt);
		LoopingUpdate_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingIntermediate_parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Intermediate(i, dt);
			}
		}, ap);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInnerInteraction2nd_parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				InnerInteraction2nd(i, dt);
			}
		}, ap);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingContactInteraction2nd_parallel_exec(Real dt)
	{
		for (size_t k = 0; k < interacting_bodies_.size(); ++k)
		{

			parallel_for(blocked_range<size_t>(0, (*indexes_interacting_particles_[k]).size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.begin(); l < r.end(); ++l) {
					size_t particle_index_i = (*indexes_interacting_particles_[k])[l];
					ContactInteraction2nd(particle_index_i, k, dt);
				}
			}, ap);
		}
	}	
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplex2Levels<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::parallel_exec(Real dt)
	{
		SetupDynamics(dt);
		LoopingInitialization_parallel_exec(dt);
		LoopingInnerInteraction_parallel_exec(dt);
		LoopingContactInteraction_parallel_exec(dt);
		LoopingIntermediate_parallel_exec(dt);
		LoopingInnerInteraction2nd_parallel_exec(dt);
		LoopingContactInteraction2nd_parallel_exec(dt);
		LoopingUpdate_parallel_exec(dt);
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInnerInteraction_exec(Real dt)
	{
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
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class InteractingBodytype, class InteractingParticlesType>
	void ParticleDynamicsComplexSplitting<BodyType, ParticlesType, InteractingBodytype, InteractingParticlesType>
		::LoopingInnerInteraction_parallel_exec(Real dt)
	{
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
	}
	//===============================================================//
}
