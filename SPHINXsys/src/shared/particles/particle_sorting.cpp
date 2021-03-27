/**
 * @file 	particle_sorting.cpp
 * @author	Xiangyu Hu
 */

#include "particle_sorting.h"

#include "base_body.h"
#include "base_particles.h"
#include "mesh_cell_linked_list.h"

namespace SPH {
	//=================================================================================================//
	SwapParticleData::SwapParticleData(BaseParticles* base_particles) :
		sequence_(base_particles->sequence_),
		unsorted_id_(base_particles->unsorted_id_),
		sortable_matrices_(base_particles->sortable_matrices_),
		sortable_vectors_(base_particles->sortable_vectors_),
		sortable_scalars_(base_particles->sortable_scalars_) {}
	//=================================================================================================//
	void SwapParticleData::operator () (size_t* a, size_t* b)
	{
		std::swap(*a, *b);

		size_t index_a = a - sequence_.data();
		size_t index_b = b - sequence_.data();
		std::swap(unsorted_id_[index_a], unsorted_id_[index_b]);
		for (size_t i = 0; i != sortable_matrices_.size(); ++i) {
			StdLargeVec<Matd>& matrices = *(sortable_matrices_[i]);
			std::swap(matrices[index_a], matrices[index_b]);
		}
		for (size_t i = 0; i != sortable_vectors_.size(); ++i) {
			StdLargeVec<Vecd>& vectors = *(sortable_vectors_[i]);
			std::swap(vectors[index_a], vectors[index_b]);
		}
		for (size_t i = 0; i != sortable_scalars_.size(); ++i) {
			StdLargeVec<Real>& scalars = *(sortable_scalars_[i]);
			std::swap(scalars[index_a], scalars[index_b]);
		}
	}	
	//=================================================================================================//
	ParticleSorting::ParticleSorting(RealBody* real_body) :
		base_particles_(NULL), swap_particle_data_(NULL), compare_(), 
		quick_sort_particle_range_(NULL), quick_sort_particle_body_() {}
	//=================================================================================================//
	void ParticleSorting::sortingParticleData(size_t* begin, size_t size)
	{
		quick_sort_particle_range_->begin_ = begin;
		quick_sort_particle_range_->size_ = size;
		parallel_for(*quick_sort_particle_range_, quick_sort_particle_body_, ap);
		updateSortedId();
	}
	//=================================================================================================//
	void ParticleSorting::updateSortedId()
	{
		StdLargeVec<size_t>& unsorted_id = base_particles_->unsorted_id_;
		StdLargeVec<size_t>& sorted_id = base_particles_->sorted_id_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		parallel_for(blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					sorted_id[unsorted_id[i]] = i;
				}
			}, ap);
	}
	//=================================================================================================//
	void ParticleSorting::assignBaseParticles(BaseParticles* base_particles)
	{
		base_particles_ = base_particles;
		swap_particle_data_ = new SwapParticleData(base_particles);
		size_t* begin = base_particles_->sequence_.data();
		quick_sort_particle_range_ = new tbb::interafce9::internal::QuickSortParticleRange<size_t*, 
			CompareParticleSequence, SwapParticleData>(begin, 0, compare_, *swap_particle_data_);
	};
	//=================================================================================================//
}
