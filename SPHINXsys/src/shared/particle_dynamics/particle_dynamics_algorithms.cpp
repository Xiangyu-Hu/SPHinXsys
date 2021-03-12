/**
* @file 	particle_dynamics_algorithms.cpp
* @brief 	This is the implementation of the template class particle dynamics algorithms
* @author	Chi ZHang and Xiangyu Hu
*/

#include "particle_dynamics_algorithms.h"

//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	void ParticleDynamicsSimple::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsSimple::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator_parallel(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void InteractionDynamics::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->exec(dt);
	}
	//=================================================================================================//
	void InteractionDynamics::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->parallel_exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator_parallel(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->parallel_exec(dt);
	}
	//=================================================================================================//
	CombinedInteractionDynamics::
		CombinedInteractionDynamics(StdVec<InteractionDynamics*> multiple_dynamics) : 
		InteractionDynamics(multiple_dynamics[0]->sph_body_),
		multiple_dynamics_(multiple_dynamics) 
	{
		for (InteractionDynamics* dynamics : multiple_dynamics) 
		{
			if (sph_body_ != dynamics->sph_body_)
			{
				std::cout << "\n Error: CombinedInteractionDynamics does not have the same source body!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
			for (size_t k = 0; k < dynamics->pre_processes_.size(); ++k)
				pre_processes_.push_back(dynamics->pre_processes_[k]);
			for (size_t k = 0; k < dynamics->post_processes_.size(); ++k)
				post_processes_.push_back(dynamics->post_processes_[k]);
		}
	}
	//=================================================================================================//
	void CombinedInteractionDynamics::setupDynamics(Real dt)
	{
		for (InteractionDynamics* dynamics : multiple_dynamics_) dynamics->setupDynamics(dt);
	}		
	//=================================================================================================//
	void CombinedInteractionDynamics::Interaction(size_t index_i, Real dt)
	{
		for (InteractionDynamics* dynamics : multiple_dynamics_) dynamics->Interaction(index_i, dt);
	}		
	//=================================================================================================//
	void InteractionDynamicsWithUpdate::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->exec(dt);
		ParticleIterator(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void InteractionDynamicsWithUpdate::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->parallel_exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator_parallel(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->parallel_exec(dt);
		ParticleIterator_parallel(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamics1Level::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator(total_real_particles, functor_initialization_, dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->exec(dt);
		ParticleIterator(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->exec(dt);
		ParticleIterator(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamics1Level::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator_parallel(total_real_particles, functor_initialization_, dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->parallel_exec(dt);
		ParticleIterator_parallel(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->parallel_exec(dt);
		ParticleIterator_parallel(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void InteractionDynamicsSplitting::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->exec(dt);
		ParticleIteratorSplittingSweep(split_cell_lists_, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->exec(dt);
	}
	//=================================================================================================//
	void InteractionDynamicsSplitting::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->parallel_exec(dt);
		ParticleIteratorSplittingSweep_parallel(split_cell_lists_, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->parallel_exec(dt);
	}
	//=============================================================================================//
}
//=================================================================================================//