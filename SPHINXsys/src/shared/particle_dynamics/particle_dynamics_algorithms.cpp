/**
* @file 	particle_dynamics_algorithms.cpp
* @brief 	This is the implementation of the template class particle dynamics algorithms
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#include "particle_dynamics_algorithms.h"

//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	void ParticleDynamicsSimple::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamicsSimple::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator_parallel(number_of_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void InteractionDynamics::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator(number_of_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->exec(dt);
	}
	//=================================================================================================//
	void InteractionDynamics::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->parallel_exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator_parallel(number_of_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->parallel_exec(dt);
	}
	//=================================================================================================//
	void InteractionDynamicsWithUpdate::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator(number_of_particles, functor_interaction_, dt);
		ParticleIterator(number_of_particles, functor_update_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->exec(dt);
	}
	//=================================================================================================//
	void InteractionDynamicsWithUpdate::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->parallel_exec(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator_parallel(number_of_particles, functor_interaction_, dt);
		ParticleIterator_parallel(number_of_particles, functor_update_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->parallel_exec(dt);
	}
	//=================================================================================================//
	void ParticleDynamics1Level::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator(number_of_particles, functor_initialization_, dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->exec(dt);
		ParticleIterator(number_of_particles, functor_interaction_, dt);
		ParticleIterator(number_of_particles, functor_update_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->exec(dt);
	}
	//=================================================================================================//
	void ParticleDynamics1Level::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t number_of_particles = sph_body_->number_of_particles_;
		ParticleIterator_parallel(number_of_particles, functor_initialization_, dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k) pre_processes_[k]->parallel_exec(dt);
		ParticleIterator_parallel(number_of_particles, functor_interaction_, dt);
		ParticleIterator_parallel(number_of_particles, functor_update_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k) post_processes_[k]->parallel_exec(dt);
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