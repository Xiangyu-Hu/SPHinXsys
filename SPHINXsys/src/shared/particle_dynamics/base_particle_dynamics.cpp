/**
 * @file 	base_particle_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */
#include "base_particle_dynamics.h"
#include "cell_linked_list.h"
//=============================================================================================//
namespace SPH
{
	Real GlobalStaticVariables::physical_time_ = 0.0;
	//=============================================================================================//
	void ParticleIterator(size_t total_real_particles, ParticleFunctor &particle_functor, Real dt)
	{
		for (size_t i = 0; i < total_real_particles; ++i)
			particle_functor(i, dt);
	}
	//=============================================================================================//
	void ParticleIterator_parallel(size_t total_real_particles, ParticleFunctor &particle_functor, Real dt)
	{
		parallel_for(blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				particle_functor(i, dt);
			}
		}, ap);
	}
	//=================================================================================================//
	void ParticleIteratorSplittingSweep(SplitCellLists& split_cell_lists,
		ParticleFunctor& particle_functor, Real dt)
	{
		Real dt2 = dt * 0.5;
		//forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			ConcurrentCellLists& cell_lists = split_cell_lists[k];
			for (size_t l = 0; l != cell_lists.size(); ++l)
			{
				IndexVector& particle_indexes
					= cell_lists[l]->real_particle_indexes_;
				for (size_t i = 0; i != particle_indexes.size(); ++i)
				{
					particle_functor(particle_indexes[i], dt2);
				}
			}
		}

		//backward sweeping
		for (size_t k = split_cell_lists.size(); k != 0; --k) {
			ConcurrentCellLists& cell_lists = split_cell_lists[k - 1];
			for (size_t l = 0; l != cell_lists.size(); ++l)
			{
				IndexVector& particle_indexes
					= cell_lists[l]->real_particle_indexes_;
				for (size_t i = particle_indexes.size(); i != 0; --i)
				{
					particle_functor(particle_indexes[i - 1], dt2);
				}
			}
		}
	}
	//=================================================================================================//
	void ParticleIteratorSplittingSweep_parallel(SplitCellLists& split_cell_lists,
		ParticleFunctor &particle_functor, Real dt)
	{
		Real dt2 = dt * 0.5;
		//forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			ConcurrentCellLists& cell_lists = split_cell_lists[k];
			parallel_for(blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t>& r) {
					for (size_t l = r.begin(); l < r.end(); ++l) {
						IndexVector& particle_indexes
							= cell_lists[l]->real_particle_indexes_;
						for (size_t i = 0; i < particle_indexes.size(); ++i)
						{
							particle_functor(particle_indexes[i], dt2);
						}
					}
				}, ap);
		}

		//backward sweeping
		for (size_t k = split_cell_lists.size(); k != 0; --k) {
			ConcurrentCellLists& cell_lists = split_cell_lists[k - 1];
			parallel_for(blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.begin(); l < r.end(); ++l) {
					IndexVector& particle_indexes
						= cell_lists[l]->real_particle_indexes_;
					for (size_t i = particle_indexes.size(); i != 0; --i)
					{
						particle_functor(particle_indexes[i - 1], dt2);
					}
				}
			}, ap);
		}
	}
	//=============================================================================================//
}
//=============================================================================================//