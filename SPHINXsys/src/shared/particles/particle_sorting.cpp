/**
 * @file 	particle_sorting.cpp
 * @author	Xiangyu Hu
 */

#include "particle_sorting.h"

#include "base_body.h"
#include "base_particles.h"
#include "base_particle_dynamics.h"
#include "cell_linked_list.h"

namespace SPH
{
	//=================================================================================================//
	SwapSortableParticleData::SwapSortableParticleData(BaseParticles *base_particles)
		: sequence_(base_particles->sequence_),
		  unsorted_id_(base_particles->unsorted_id_),
		  sortable_data_(base_particles->sortable_data_) {}
	//=================================================================================================//
	void SwapSortableParticleData::operator()(size_t *a, size_t *b)
	{
		std::swap(*a, *b);

		size_t index_a = a - sequence_.data();
		size_t index_b = b - sequence_.data();
		std::swap(unsorted_id_[index_a], unsorted_id_[index_b]);
		swap_particle_data_value_(sortable_data_, index_a, index_b);
	}
	//=================================================================================================//
	ParticleSorting::ParticleSorting(RealBody *real_body)
		: base_particles_(nullptr), swap_sortable_particle_data_(nullptr), compare_(),
		  quick_sort_particle_range_(nullptr), quick_sort_particle_body_() {}
	//=================================================================================================//
	void ParticleSorting::sortingParticleData(size_t *begin, size_t size)
	{
		quick_sort_particle_range_->begin_ = begin;
		quick_sort_particle_range_->size_ = size;
		parallel_for(*quick_sort_particle_range_, quick_sort_particle_body_, ap);
		updateSortedId();
	}
	//=================================================================================================//
	void ParticleSorting::updateSortedId()
	{
		const StdLargeVec<size_t> &unsorted_id = base_particles_->unsorted_id_;
		StdLargeVec<size_t> &sorted_id = base_particles_->sorted_id_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		parallel_for(
			blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					sorted_id[unsorted_id[i]] = i;
				}
			},
			ap);
	}
	//=================================================================================================//
	void ParticleSorting::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
		swap_sortable_particle_data_ = swap_particle_ptr_keeper_.createPtr<SwapSortableParticleData>(base_particles);
		size_t *begin = base_particles_->sequence_.data();
		quick_sort_particle_range_ =
			quick_sort_particle_range_ptr_keeper_.createPtr<
				tbb::interface9::QuickSortParticleRange<
					size_t *, CompareParticleSequence, SwapSortableParticleData>>(begin, 0, compare_, *swap_sortable_particle_data_);
	};
	//=================================================================================================//
}
