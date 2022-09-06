/**
 * @file 	particle_dynamics_algorithms.cpp
 * @brief 	This is the implementation of the template class particle dynamics algorithms
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "particle_dynamics_algorithms.h"
#include "base_body.h"
#include "base_particles.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	void OldInteractionDynamics::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->exec(dt);
	}
	//=================================================================================================//
	void OldInteractionDynamics::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->parallel_exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator_parallel(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->parallel_exec(dt);
	}
	//=================================================================================================//
	void OldInteractionDynamicsWithUpdate::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->exec(dt);
		ParticleIterator(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void OldInteractionDynamicsWithUpdate::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->parallel_exec(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator_parallel(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->parallel_exec(dt);
		ParticleIterator_parallel(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamics1Level::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator(total_real_particles, functor_initialization_, dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->exec(dt);
		ParticleIterator(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->exec(dt);
		ParticleIterator(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	void ParticleDynamics1Level::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		size_t total_real_particles = base_particles_->total_real_particles_;
		ParticleIterator_parallel(total_real_particles, functor_initialization_, dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->parallel_exec(dt);
		ParticleIterator_parallel(total_real_particles, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->parallel_exec(dt);
		ParticleIterator_parallel(total_real_particles, functor_update_, dt);
	}
	//=================================================================================================//
	InteractionDynamicsSplitting::InteractionDynamicsSplitting(SPHBody &sph_body)
		: OldInteractionDynamics(sph_body),
		  real_body_(DynamicCast<RealBody>(this, sph_body)),
		  split_cell_lists_(real_body_.getSplitCellLists())
	{
		real_body_.setUseSplitCellLists();
	};
	//=================================================================================================//
	void InteractionDynamicsSplitting::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->exec(dt);
		ParticleIteratorSplittingSweep(split_cell_lists_, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->exec(dt);
	}
	//=================================================================================================//
	void InteractionDynamicsSplitting::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t k = 0; k < pre_processes_.size(); ++k)
			pre_processes_[k]->parallel_exec(dt);
		ParticleIteratorSplittingSweep_parallel(split_cell_lists_, functor_interaction_, dt);
		for (size_t k = 0; k < post_processes_.size(); ++k)
			post_processes_[k]->parallel_exec(dt);
	}
	//=============================================================================================//
}
//=================================================================================================//