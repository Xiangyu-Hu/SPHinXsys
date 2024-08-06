#include "particle_sorting.h"

#include "base_body.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
SwapSortableParticleData::SwapSortableParticleData(BaseParticles &base_particles)
    : sequence_(base_particles.ParticleSequences()),
      sortable_data_(base_particles.SortableParticleData()),
      swap_particle_data_value_(sortable_data_) {}
//=================================================================================================//
void SwapSortableParticleData::operator()(size_t *a, size_t *b)
{
    std::swap(*a, *b);

    size_t index_a = a - sequence_.data();
    size_t index_b = b - sequence_.data();
    swap_particle_data_value_(index_a, index_b);
}
//=================================================================================================//
ParticleSorting::ParticleSorting(BaseParticles &base_particles)
    : base_particles_(base_particles),
      original_id_(base_particles.ParticleOriginalIds()),
      sorted_id_(base_particles.ParticleSortedIds()),
      sequence_(base_particles.ParticleSequences()),
      swap_sortable_particle_data_(base_particles), compare_(),
      quick_sort_particle_range_(sequence_.data(), 0, compare_, swap_sortable_particle_data_),
      quick_sort_particle_body_()
{
    base_particles.addVariableToSort<size_t>("OriginalID");
}
//=================================================================================================//
void ParticleSorting::sortingParticleData(size_t *begin, size_t size)
{
    quick_sort_particle_range_.begin_ = begin;
    quick_sort_particle_range_.size_ = size;
    parallel_for(quick_sort_particle_range_, quick_sort_particle_body_, ap);
    updateSortedId();
}
//=================================================================================================//
void ParticleSorting::updateSortedId()
{
    size_t total_real_particles = base_particles_.TotalRealParticles();
    parallel_for(
        IndexRange(0, total_real_particles),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                sorted_id_[original_id_[i]] = i;
            }
        },
        ap);
}
//=================================================================================================//
} // namespace SPH
