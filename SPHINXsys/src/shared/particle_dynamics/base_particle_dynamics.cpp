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
	void ContactIterator(StdVec<ListIndexVector*> indexes_interacting_particles,
		ContactFunctor &contact_functor, Real dt)
	{
		for (size_t k = 0; k < indexes_interacting_particles.size(); ++k)
			for (size_t l = 0; l < (*indexes_interacting_particles[k]).size(); ++l) {
				size_t particle_index_i = (*indexes_interacting_particles[k])[l];
				contact_functor(particle_index_i, k, dt);
			}
	}
	//===============================================================//
	void ContactIterator_parallel(StdVec<ListIndexVector*> indexes_interacting_particles,
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
	//===============================================================//
	void InnerIteratorSplitting(ByCellLists by_cell_lists_particle_indexes,
		InnerFunctor &inner_functor, Real dt)
	{
		size_t number_of_lists = (*by_cell_lists_particle_indexes).size();

		//forward sweeping
		for (size_t k = 0; k < number_of_lists; ++k) {
			StdVec<IndexVector> &lists_particle_indexes
				= by_cell_lists_particle_indexes[k];
			for (size_t l = 0; l < lists_particle_indexes.size(); ++l)
			{
				for (size_t i = 0; i < lists_particle_indexes[l].size(); ++i)
				{
					inner_functor(lists_particle_indexes[l][i], dt);
				}
			}
		}

		//backward sweeping
		for (size_t k = number_of_lists - 1; k >= 0; --k) {
			StdVec<IndexVector> &lists_particle_indexes
				= by_cell_lists_particle_indexes[k];
			for (size_t l = lists_particle_indexes.size() - 1; l >= 0; --l)
			{
				for (size_t i = lists_particle_indexes[l].size() - 1; i >= 0; --i)
				{
					inner_functor(lists_particle_indexes[l][i], dt);
				}
			}
		}
	}
	//===============================================================//
	void InnerIteratorSplitting_parallel(ByCellLists by_cell_lists_particle_indexes,
		InnerFunctor &inner_functor, Real dt)
	{
		size_t number_of_lists = (*by_cell_lists_particle_indexes).size();
		//forward sweeping
		for (size_t k = 0; k < number_of_lists; ++k) {
			StdVec<IndexVector> &lists_particle_indexes
				= by_cell_lists_particle_indexes[k];
			parallel_for(blocked_range<size_t>(0, lists_particle_indexes.size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.begin(); l < r.end(); ++l) {
					for (size_t i = 0; i < lists_particle_indexes[l].size(); ++i)
					{
						inner_functor(lists_particle_indexes[l][i], dt);
					}
				}
			}, ap);
		}

		//backward sweeping
		for (size_t k = number_of_lists; k >= 1; --k) {
			StdVec<IndexVector> &lists_particle_indexes
				= by_cell_lists_particle_indexes[k - 1];
			parallel_for(blocked_range<size_t>(0, lists_particle_indexes.size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t l = r.end(); l >= r.begin() + 1; --l) {
					IndexVector &particle_indexes = lists_particle_indexes[l - 1];
					for (size_t i = particle_indexes.size(); i >= 1; --i)
					{
						inner_functor(particle_indexes[i - 1], dt);
					}
				}
			}, ap);
		}
	}
}