/**
 * @file 	base_particle_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH
{
	Real GlobalStaticVariables::physical_time_ = 0.0;

	//===============================================================//
	void InnerIterator(size_t number_of_particles, InnerFunctor &inner_functor, Real dt)
	{
		for (size_t i = 0; i < number_of_particles; ++i)
			inner_functor(i, dt);
	}
	//===============================================================//
	void InnerIterator_parallel(size_t number_of_particles, InnerFunctor &inner_functor, Real dt)
	{
		parallel_for(blocked_range<size_t>(0, number_of_particles),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				inner_functor(i, dt);
			}
		}, ap);
	}
	//===============================================================//
	void ContactIterator(InteractingParticles& indexes_interacting_particles,
		ContactFunctor &contact_functor, Real dt)
	{
		for (size_t k = 0; k < indexes_interacting_particles.size(); ++k)
			for (size_t l = 0; l < (*indexes_interacting_particles[k]).size(); ++l) {
				size_t particle_index_i = (*indexes_interacting_particles[k])[l];
				contact_functor(particle_index_i, k, dt);
			}
	}
	//===============================================================//
	void ContactIterator_parallel(InteractingParticles& indexes_interacting_particles,
		ContactFunctor &contact_functor, Real dt)
	{
		for (size_t k = 0; k < indexes_interacting_particles.size(); ++k)
			parallel_for(blocked_range<size_t>(0, (*indexes_interacting_particles[k]).size()),
				[&](const blocked_range<size_t>& r) {
			for (size_t l = r.begin(); l < r.end(); ++l) {
				size_t particle_index_i = (*indexes_interacting_particles[k])[l];
				contact_functor(particle_index_i, k, dt);
			}
		}, ap);
	}
	//=================================================================================================//
	void CellListIteratorSplitting(SplitCellLists& split_cell_lists,
		CellListFunctor& cell_list_functor, Real dt)
	{
		//forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			for (size_t l = 0; l != cell_lists.size(); ++l)
				cell_list_functor(cell_lists[l], dt);
		}

		//backward sweeping
		for (size_t k = split_cell_lists.size(); k != 0; --k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k - 1];
			for (size_t l = cell_lists.size(); l != 0; --l)
				cell_list_functor(cell_lists[l - 1], dt);
		}
	}
	//=================================================================================================//
	void CellListIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		CellListFunctor& cell_list_functor, Real dt)
	{
		//forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			parallel_for(blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t>& r) {
					for (size_t l = r.begin(); l < r.end(); ++l)
						cell_list_functor(cell_lists[l], dt);
				}, ap);
		}
	
		//backward sweeping
		for (size_t k = split_cell_lists.size(); k >= 1; --k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k - 1];
			parallel_for(blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t>& r) {
					for (size_t l = r.end(); l >= r.begin() + 1; --l) {
						cell_list_functor(cell_lists[l -1], dt);
					}
				}, ap);
		}
	}
	//===============================================================//
	void InnerIteratorSplitting(SplitCellLists& split_cell_lists,
		InnerFunctor &inner_functor, Real dt)
	{
		//forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			for (size_t l = 0; l != cell_lists.size(); ++l)
			{
				ConcurrentListDataVector& particle_data_lists
					= cell_lists[l]->particle_data_lists_;
				for (size_t i = 0; i != particle_data_lists.size(); ++i)
				{
					inner_functor(particle_data_lists[i].first, dt);
				}
			}
		}

		//backward sweeping
		for (size_t k = split_cell_lists.size() - 1; k >= 0; --k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			for (size_t l = cell_lists.size() - 1; l >= 0; --l)
			{
				ConcurrentListDataVector& particle_data_lists
					= cell_lists[l]->particle_data_lists_;
				for (size_t i = particle_data_lists.size() - 1; i >= 0; --i)
				{
					inner_functor(particle_data_lists[i].first, dt);
				}
			}
		}
	}
	//===============================================================//
	void InnerIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		InnerFunctor &inner_functor, Real dt)
	{
		//forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k];
			parallel_for(blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t>& r) {
					for (size_t l = r.begin(); l < r.end(); ++l) {
						ConcurrentListDataVector& particle_data_lists
							= cell_lists[l]->particle_data_lists_;
						for (size_t i = 0; i < particle_data_lists.size(); ++i)
						{
							inner_functor(particle_data_lists[i].first, dt);
						}
					}
				}, ap);
		}

		//backward sweeping
		for (size_t k = split_cell_lists.size(); k >= 1; --k) {
			StdLargeVec<CellList*>& cell_lists = split_cell_lists[k - 1];
			parallel_for(blocked_range<size_t>(0, cell_lists.size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.end(); l >= r.begin() + 1; --l) {
					ConcurrentListDataVector& particle_data_lists
						= cell_lists[l-1]->particle_data_lists_;
					for (size_t i = particle_data_lists.size(); i >= 1; --i)
					{
						inner_functor(particle_data_lists[i - 1].first, dt);
					}
				}
			}, ap);
		}
	}
	//===============================================================//
}