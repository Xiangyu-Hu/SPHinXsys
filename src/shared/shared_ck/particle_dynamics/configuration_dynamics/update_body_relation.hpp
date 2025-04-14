#ifndef UPDATE_BODY_RELATION_HPP
#define UPDATE_BODY_RELATION_HPP

#include "update_body_relation.h"

#include "cell_linked_list.hpp"
#include "particle_iterators_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    UpdateRelation(Relation<Inner<Parameters...>> &inner_relation)
    : Interaction<Inner<Parameters...>>(inner_relation),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{}),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, this->real_body_->getCellLinkedList())),
      kernel_implementation_(*this) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
template <class EncloserType>
UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : Interaction<Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      neighbor_search_(encloser.cell_linked_list_.createNeighborSearch(ex_policy)),
      grid_spacing_squared_(
          pow(encloser.cell_linked_list_.getMesh().GridSpacing(), 2)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    InteractKernel::incrementNeighborSize(UnsignedInt index_i)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, this->source_pos_,
        [&](size_t index_j)
        {
            if (index_i != index_j)
            {
                if ((this->source_pos_[index_i] - this->target_pos_[index_j])
                        .squaredNorm() < grid_spacing_squared_)
                    neighbor_count++;
            }
        });
    this->neighbor_index_[index_i] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    InteractKernel::updateNeighborList(UnsignedInt index_i)
{
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        index_i, this->source_pos_,
        [&](size_t index_j)
        {
            if (index_i != index_j)
            {
                if ((this->source_pos_[index_i] - this->target_pos_[index_j])
                        .squaredNorm() < grid_spacing_squared_)
                {
                    this->neighbor_index_[this->particle_offset_[index_i] + neighbor_count] = index_j;
                    neighbor_count++;
                }
            }
        });
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();
    InteractKernel *computing_kernel = kernel_implementation_.getComputingKernel();
    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementNeighborSize(i); });

    auto *dv_neighbor_index = this->inner_relation_.getNeighborIndex();
    auto *dv_particle_offset = this->inner_relation_.getParticleOffset();
    UnsignedInt *neighbor_index = dv_neighbor_index->DelegatedData(ex_policy_);
    UnsignedInt *particle_offset = dv_particle_offset->DelegatedData(ex_policy_);
    UnsignedInt current_offset_list_size = total_real_particles + 1;
    UnsignedInt current_neighbor_index_size =
        exclusive_scan(ex_policy_, neighbor_index, particle_offset, current_offset_list_size,
                       typename PlusUnsignedInt<ExecutionPolicy>::type());

    if (current_neighbor_index_size > dv_neighbor_index->getDataSize())
    {
        dv_neighbor_index->reallocateData(ex_policy_, current_neighbor_index_size);
        this->inner_relation_.resetComputingKernelUpdated();
        kernel_implementation_.overwriteComputingKernel();
    }

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborList(i); });
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    UpdateRelation(RelationType &contact_relation)
    : Interaction<Contact<Parameters...>>(contact_relation),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{})
{
    for (auto &contact_body : this->contact_bodies_)
    {
        contact_cell_linked_list_.push_back(
            DynamicCast<CellLinkedList>(this, &contact_body->getCellLinkedList()));
        contact_kernel_implementation_.push_back(
            contact_kernel_implementation_ptrs_.template createPtr<KernelImplementation>(*this));
    }
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
template <class EncloserType>
UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : Interaction<Contact<Parameters...>>::InteractKernel(ex_policy, encloser, contact_index),
      masked_criterion_(
          ex_policy, encloser.contact_relation_.getContactIdentifier(contact_index), *this),
      neighbor_search_(
          encloser.contact_cell_linked_list_[contact_index]->createNeighborSearch(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::incrementNeighborSize(UnsignedInt source_index)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        source_index, this->source_pos_,
        [&](size_t target_index)
        {
            if (masked_criterion_(target_index, source_index))
                neighbor_count++;
        });
    this->neighbor_index_[source_index] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::updateNeighborList(UnsignedInt source_index)
{
    UnsignedInt neighbor_count = 0;
    neighbor_search_.forEachSearch(
        source_index, this->source_pos_,
        [&](size_t target_index)
        {
            if (masked_criterion_(target_index, source_index))
            {
                this->neighbor_index_[this->particle_offset_[source_index] + neighbor_count] = target_index;
                neighbor_count++;
            }
        });
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();

    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        InteractKernel *computing_kernel = contact_kernel_implementation_[k]->getComputingKernel(k);
        particle_for(ex_policy_,
                     IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { computing_kernel->incrementNeighborSize(i); });

        auto *dv_neighbor_index = this->contact_relation_.getNeighborIndex(k);
        auto *dv_particle_offset = this->contact_relation_.getParticleOffset(k);
        UnsignedInt *neighbor_index = dv_neighbor_index->DelegatedData(ex_policy_);
        UnsignedInt *particle_offset = dv_particle_offset->DelegatedData(ex_policy_);
        UnsignedInt current_offset_list_size = total_real_particles + 1;
        UnsignedInt current_neighbor_index_size =
            exclusive_scan(ex_policy_, neighbor_index, particle_offset, current_offset_list_size,
                           typename PlusUnsignedInt<ExecutionPolicy>::type());

        if (current_neighbor_index_size > dv_neighbor_index->getDataSize())
        {
            dv_neighbor_index->reallocateData(ex_policy_, current_neighbor_index_size);
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
template <class ExecutionPolicy, class FirstRelation, class... Others>
template <class FirstParameterSet, typename... OtherParameterSets>
UpdateRelation<ExecutionPolicy, FirstRelation, Others...>::UpdateRelation(
    FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets)
    : UpdateRelation<ExecutionPolicy, FirstRelation>(first_parameter_set),
      other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...) {}
//=================================================================================================//
template <class ExecutionPolicy, class FirstRelation, class... Others>
void UpdateRelation<ExecutionPolicy, FirstRelation, Others...>::exec(Real dt)
{
    UpdateRelation<ExecutionPolicy, FirstRelation>::exec(dt);
    other_interactions_.exec(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_BODY_RELATION_HPP
