#ifndef UPDATE_BODY_RELATION_HPP
#define UPDATE_BODY_RELATION_HPP

#include "update_body_relation.h"

#include "cell_linked_list.hpp"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
BodyRelationUpdate<Inner<>>::ComputingKernel<ExecutionPolicy>::ComputingKernel(
    const ExecutionPolicy &ex_policy, BodyRelationUpdate<Inner<>> &encloser)
    : LocalInteractionDynamics<Inner<>>::ComputingKernel<ExecutionPolicy>(ex_policy, encloser),
      neighbor_search_(encloser.cell_linked_list_.createNeighborSearch(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
void BodyRelationUpdate<Inner<>>::ComputingKernel<ExecutionPolicy>::
    incrementNeighborSize(UnsignedInt index_i)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, this->source_pos_,
        [&](size_t index_j)
        {
            if (index_i != index_j)
            {
                neighbor_count++;
            }
        });
    this->neighbor_index_[index_i] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy>
void BodyRelationUpdate<Inner<>>::ComputingKernel<ExecutionPolicy>::
    updateNeighborList(UnsignedInt index_i)
{
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, this->source_pos_,
        [&](size_t index_j)
        {
            if (index_i != index_j)
            {
                neighbor_count++;
            }
            this->neighbor_index_[this->particle_offset_[index_i] + neighbor_count] = index_j;
        });
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
template <typename... Args>
UpdateRelation<ExecutionPolicy, BodyRelationUpdate<Inner<Parameters...>>>::UpdateRelation(Args &&...args)
    : BodyRelationUpdate<Inner<Parameters...>>(std::forward<Args>(args)...),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{}), kernel_implementation_(*this) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, BodyRelationUpdate<Inner<Parameters...>>>::exec(Real dt)
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
        exclusive_scan(ex_policy_, neighbor_index, particle_offset,
                       this->particle_offset_list_size_,
                       typename PlusUnsignedInt<ExecutionPolicy>::type());

    if (current_neighbor_index_size > this->dv_neighbor_index_->getDataFieldSize())
    {
        this->dv_neighbor_index_->reallocateDataField(ex_policy_, current_neighbor_index_size);
        this->inner_relation_.resetComputingKernelUpdated();
        kernel_implementation_.overwriteComputingKernel();
    }

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborList(i); });
}
//=================================================================================================//
template <class ExecutionPolicy>
BodyRelationUpdate<Contact<>>::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    BodyRelationUpdate<Contact<>> &encloser,
                    UnsignedInt contact_index)
    : LocalInteractionDynamics<Contact<>>::ComputingKernel<ExecutionPolicy>(
          ex_policy, encloser, contact_index),
      neighbor_search_(encloser.contact_cell_linked_list_[contact_index]
                           ->createNeighborSearch(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
void BodyRelationUpdate<Contact<>>::ComputingKernel<ExecutionPolicy>::
    incrementNeighborSize(UnsignedInt index_i)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, this->source_pos_,
        [&](size_t index_j)
        {
            neighbor_count++;
        });
    this->neighbor_index_[index_i] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy>
void BodyRelationUpdate<Contact<>>::ComputingKernel<ExecutionPolicy>::
    updateNeighborList(UnsignedInt index_i)
{
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, this->source_pos_,
        [&](size_t index_j)
        {
            neighbor_count++;
            this->neighbor_index_[this->particle_offset_[index_i] + neighbor_count] = index_j;
        });
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
template <typename... Args>
UpdateRelation<ExecutionPolicy, BodyRelationUpdate<Contact<Parameters...>>>::UpdateRelation(Args &&...args)
    : BodyRelationUpdate<Contact<Parameters...>>(std::forward<Args>(args)...),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{})
{
    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        contact_kernel_implementation_.push_back(
            contact_kernel_implementation_ptrs_.template createPtr<KernelImplementation>(*this));
    }
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, BodyRelationUpdate<Contact<Parameters...>>>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();

    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        ComputingKernel *computing_kernel = contact_kernel_implementation_[k]->getComputingKernel(k);
        particle_for(ex_policy_,
                     IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { computing_kernel->incrementNeighborSize(i); });

        UnsignedInt *neighbor_index = this->dv_contact_neighbor_index_[k]->DelegatedDataField(ex_policy_);
        UnsignedInt *particle_offset = this->dv_contact_particle_offset_[k]->DelegatedDataField(ex_policy_);
        UnsignedInt current_neighbor_index_size =
            exclusive_scan(ex_policy_, neighbor_index, particle_offset,
                           this->particle_offset_list_size_,
                           typename PlusUnsignedInt<ExecutionPolicy>::type());

        if (current_neighbor_index_size > this->dv_contact_neighbor_index_[k]->getDataFieldSize())
        {
            this->dv_contact_neighbor_index_[k]->reallocateDataField(ex_policy_, current_neighbor_index_size);
            this->contact_relation_.resetComputingKernelUpdated(k);
            contact_kernel_implementation_[k]->overwriteComputingKernel(k);
        }

        particle_for(ex_policy_,
                     IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { computing_kernel->updateNeighborList(i); });
    }
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_BODY_RELATION_HPP
