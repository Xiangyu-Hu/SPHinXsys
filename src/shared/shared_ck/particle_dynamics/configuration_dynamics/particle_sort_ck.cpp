#include "particle_sort_ck.hpp"

namespace SPH
{
//=================================================================================================//
UpdateSortableVariables::UpdateSortableVariables(BaseParticles *particles)
    : initialize_temp_variables_()
{
    initialize_temp_variables_(temp_variables_, particles->ParticlesBound());
}
//=================================================================================================//
QuickSort::SwapParticleIndex::SwapParticleIndex(UnsignedInt *sequence, UnsignedInt *index_permutation)
    : sequence_(sequence), index_permutation_(index_permutation) {}
//=================================================================================================//
void QuickSort::SwapParticleIndex::operator()(UnsignedInt *a, UnsignedInt *b)
{
    std::swap(*a, *b);

    UnsignedInt index_a = a - sequence_;
    UnsignedInt index_b = b - sequence_;
    std::swap(index_permutation_[index_a], index_permutation_[index_b]);
}
//=================================================================================================//
void QuickSort::sort(const ParallelPolicy &ex_policy, BaseParticles *particles)
{
    quick_sort_particle_range_.begin_ = sequence_;
    quick_sort_particle_range_.size_ = particles->TotalRealParticles();
    tbb::parallel_for(quick_sort_particle_range_, quick_sort_particle_body_);
}
//=================================================================================================//
} // namespace SPH