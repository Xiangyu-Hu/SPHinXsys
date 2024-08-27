#ifndef UPDATE_BODY_RELATION_HPP
#define UPDATE_BODY_RELATION_HPP

#include "update_body_relation.h"

#include "cell_linked_list.hpp"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
BodyRelationUpdate<Inner<>>::ComputingKernel<ExecutionPolicy>::ComputingKernel(
    const ExecutionPolicy &ex_policy, BodyRelationUpdate<Inner<>> &update_inner_relation)
    : neighbor_search_(update_inner_relation.cell_linked_list_.createNeighborSearch(ex_policy)),
      pos_(update_inner_relation.dv_pos_->DelegatedDataField(ex_policy)),
      neighbor_index_(update_inner_relation.dv_neighbor_index_->DelegatedDataField(ex_policy)),
      particle_offset_(update_inner_relation.dv_particle_offset_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
void BodyRelationUpdate<Inner<>>::ComputingKernel<ExecutionPolicy>::
    incrementNeighborSize(UnsignedInt index_i)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, pos_,
        [&](size_t index_j)
        {
            if (index_i != index_j)
            {
                neighbor_count++;
            }
        });
    neighbor_index_[index_i] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy>
void BodyRelationUpdate<Inner<>>::ComputingKernel<ExecutionPolicy>::
    updateNeighborList(UnsignedInt index_i)
{
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, pos_,
        [&](size_t index_j)
        {
            if (index_i != index_j)
            {
                neighbor_count++;
            }
            neighbor_index_[particle_offset_[index_i] + neighbor_count] = index_j;
        });
}
//=================================================================================================//
template <class BodyRelationUpdateType, class ExecutionPolicy>
template <typename... Args>
UpdateRelation<BodyRelationUpdateType, ExecutionPolicy>::UpdateRelation(Args &&...args)
    : BodyRelationUpdateType(std::forward<Args>(args)...),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{}), kernel_implementation_(*this) {}
//=================================================================================================//
template <class BodyRelationUpdateType, class ExecutionPolicy>
void UpdateRelation<BodyRelationUpdateType, ExecutionPolicy>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();
    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementNeighborSize(i); });

    UnsignedInt *neighbor_index = this->dv_neighbor_index_->DelegatedDataField(ex_policy_);
    UnsignedInt *particle_offset = this->dv_particle_offset_->DelegatedDataField(ex_policy_);
    UnsignedInt current_neighbor_index_size =
        exclusive_scan(ex_policy_, neighbor_index,
                       neighbor_index + this->particle_offset_list_size_,
                       particle_offset,
                       typename PlusUnsignedInt<ExecutionPolicy>::type());

    if (current_neighbor_index_size > this->dv_neighbor_index_->getDataFieldSize())
    {
        this->dv_neighbor_index_->reallocateDataField(current_neighbor_index_size);
        kernel_implementation_.overwriteComputingKernel();
    }

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborList(i); });
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_BODY_RELATION_HPP
