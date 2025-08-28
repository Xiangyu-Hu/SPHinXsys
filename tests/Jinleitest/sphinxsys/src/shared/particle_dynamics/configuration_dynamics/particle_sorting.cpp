#include "particle_sorting.h"

#include "base_body.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
SwapSortableParticleData::SwapSortableParticleData(BaseParticles *base_particles)
    : sequence_(base_particles->getVariableDataByName<UnsignedInt>("Sequence")),
      evolving_variables_data_(base_particles->EvolvingVariablesData()),
      swap_particle_data_value_() {}
//=================================================================================================//
void SwapSortableParticleData::operator()(UnsignedInt *a, UnsignedInt *b)
{
    std::swap(*a, *b);

    UnsignedInt index_a = a - sequence_;
    UnsignedInt index_b = b - sequence_;
    swap_particle_data_value_(evolving_variables_data_, index_a, index_b);
}
//=================================================================================================//
ParticleSequence::ParticleSequence(RealBody &real_body)
    : LocalDynamics(real_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      sequence_(particles_->registerDiscreteVariable<UnsignedInt>("Sequence", particles_->ParticlesBound())),
      cell_linked_list_(real_body.getCellLinkedList()) {}
//=================================================================================================//
void ParticleSequence::update(size_t index_i, Real dt)
{
    sequence_[index_i] = cell_linked_list_.computingSequence(pos_[index_i], index_i);
}
//=================================================================================================//
ParticleDataSort<ParallelPolicy>::ParticleDataSort(RealBody &real_body)
    : LocalDynamics(real_body), BaseDynamics<void>(),
      sequence_(particles_->getVariableDataByName<UnsignedInt>("Sequence")),
      swap_sortable_particle_data_(particles_), compare_(),
      quick_sort_particle_range_(sequence_, 0, compare_, swap_sortable_particle_data_),
      quick_sort_particle_body_() {}
//=================================================================================================//
void ParticleDataSort<ParallelPolicy>::exec(Real dt)
{
    quick_sort_particle_range_.begin_ = sequence_;
    quick_sort_particle_range_.size_ = particles_->TotalRealParticles();
    parallel_for(quick_sort_particle_range_, quick_sort_particle_body_, ap);
}
//=================================================================================================//
UpdateSortedID::UpdateSortedID(RealBody &real_body)
    : LocalDynamics(real_body),
      original_id_(particles_->getVariableDataByName<UnsignedInt>("OriginalID")),
      sorted_id_(particles_->getVariableDataByName<UnsignedInt>("SortedID")) {}
//=================================================================================================//
void UpdateSortedID::update(size_t index_i, Real dt)
{
    sorted_id_[original_id_[index_i]] = index_i;
}
//=================================================================================================//
} // namespace SPH
