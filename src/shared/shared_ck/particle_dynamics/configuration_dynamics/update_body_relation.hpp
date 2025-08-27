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
    UpdateRelation(Inner<Parameters...> &inner_relation)
    : BaseLocalDynamicsType(inner_relation.getDynamicsIdentifier()),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{}),
      inner_relation_(inner_relation),
      cell_linked_list_(DynamicCast<CellLinkedList>(
          this, inner_relation.getDynamicsIdentifier().getCellLinkedList())),
      kernel_implementation_(*this) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
template <class EncloserType>
UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : NeighborList(ex_policy, encloser.inner_relation_),
      source_pos_(encloser.inner_relation_.getSourcePosition()->DelegatedData(ex_policy)),
      masked_source_(ex_policy, encloser.inner_relation_.getDynamicsIdentifier()),
      masked_criterion_(
          ex_policy, encloser.inner_relation_.getDynamicsIdentifier(),
          ex_policy, encloser.inner_relation_.getNeighborhood()),
      neighbor_search_(encloser.cell_linked_list_.createNeighborSearch(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    InteractKernel::incrementNeighborSize(UnsignedInt source_index)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    if (masked_source_(source_index))
    {
        neighbor_search_.forEachSearch(
            source_index, source_pos_,
            [&](size_t target_index)
            {
                if (source_index != target_index)
                {
                    if (masked_criterion_(target_index, source_index))
                        neighbor_count++;
                }
            });
    }
    this->neighbor_index_[source_index] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    InteractKernel::updateNeighborList(UnsignedInt source_index)
{
    UnsignedInt neighbor_count = 0;
    if (masked_source_(source_index))
    {
        neighbor_search_.forEachSearch(
            source_index, source_pos_,
            [&](size_t target_index)
            {
                if (source_index != target_index)
                {
                    if (masked_criterion_(target_index, source_index))
                    {
                        this->neighbor_index_[this->particle_offset_[source_index] + neighbor_count] = target_index;
                        neighbor_count++;
                    }
                }
            });
    }
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

    this->logger_->debug("UpdateCellLinkedList: incrementNeighborSize done at {} for Relation {}.",
                         this->sph_body_.getName(), type_name<Inner<Parameters...>>());

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
        UnsignedInt old_size = dv_neighbor_index->getDataSize();
        dv_neighbor_index->reallocateData(ex_policy_, current_neighbor_index_size);
        this->inner_relation_.resetComputingKernelUpdated();
        kernel_implementation_.overwriteComputingKernel();

        this->logger_->info(
            "UpdateRelation: increase neighbor index size from {} to {} at {} .",
            old_size, dv_neighbor_index->getDataSize(), this->sph_body_.getName());
    }

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborList(i); });

    this->logger_->debug("UpdateCellLinkedList: updateNeighborList done at {} for Relation {}.",
                         this->sph_body_.getName(), type_name<Inner<Parameters...>>());
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    UpdateRelation(ContactRelationType &contact_relation)
    : BaseLocalDynamicsType(contact_relation.getSourceIdentifier()),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{}),
      contact_relation_(contact_relation)
{
    for (size_t k = 0; k != contact_relation.getContactBodies().size(); ++k)
    {
        contact_cell_linked_list_.push_back(
            DynamicCast<CellLinkedList>(
                this, &contact_relation.getContactIdentifier(k).getCellLinkedList()));
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
    : NeighborList(ex_policy, encloser.contact_relation_, contact_index),
      source_pos_(encloser.contact_relation_.getSourcePosition()->DelegatedData(ex_policy)),
      masked_source_(ex_policy, encloser.contact_relation_.getSourceIdentifier()),
      masked_criterion_(
          ex_policy, encloser.contact_relation_.getContactIdentifier(contact_index),
          ex_policy, encloser.contact_relation_.getNeighborhood(contact_index)),
      neighbor_search_(
          encloser.contact_cell_linked_list_[contact_index]->createNeighborSearch(ex_policy)),
      search_depth_(ex_policy, encloser.contact_relation_.getNeighborhood(contact_index)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::incrementNeighborSize(UnsignedInt source_index)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    if (masked_source_(source_index))
    {
        neighbor_search_.forEachSearch(
            source_index, source_pos_,
            [&](size_t target_index)
            {
                if (masked_criterion_(target_index, source_index))
                    neighbor_count++;
            },
            search_depth_(source_index));
    }
    this->neighbor_index_[source_index] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::updateNeighborList(UnsignedInt source_index)
{
    UnsignedInt neighbor_count = 0;
    if (masked_source_(source_index))
    {
        neighbor_search_.forEachSearch(
            source_index, source_pos_,
            [&](size_t target_index)
            {
                if (masked_criterion_(target_index, source_index))
                {
                    this->neighbor_index_[this->particle_offset_[source_index] + neighbor_count] = target_index;
                    neighbor_count++;
                }
            },
            search_depth_(source_index));
    }
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();

    for (size_t k = 0; k != contact_relation_.getContactBodies().size(); ++k)
    {
        InteractKernel *computing_kernel = contact_kernel_implementation_[k]->getComputingKernel(k);
        particle_for(ex_policy_,
                     IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { computing_kernel->incrementNeighborSize(i); });

        this->logger_->debug("UpdateRelation: incrementNeighborSize done at {} for Relation {} to {}.",
                             this->sph_body_.getName(), type_name<Contact<Parameters...>>(),
                             contact_relation_.getContactIdentifier(k).getName());

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
            UnsignedInt old_size = dv_neighbor_index->getDataSize();
            dv_neighbor_index->reallocateData(ex_policy_, current_neighbor_index_size);
            this->contact_relation_.resetComputingKernelUpdated(k);
            contact_kernel_implementation_[k]->overwriteComputingKernel(k);

            this->logger_->info(
                "UpdateRelation: increase neighbor index size from {} to {} at .",
                old_size, dv_neighbor_index->getDataSize(), this->sph_body_.getName());
        }

        particle_for(ex_policy_,
                     IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { computing_kernel->updateNeighborList(i); });

        this->logger_->debug(
            "UpdateRelation: updateNeighborList done at {} for Relation {} to {}.",
            this->sph_body_.getName(), type_name<Contact<Parameters...>>(),
            contact_relation_.getContactIdentifier(k).getName());
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class FirstRelation, class... OtherRelations>
UpdateRelation<ExecutionPolicy, FirstRelation, OtherRelations...>::UpdateRelation(
    FirstRelation &first_relation, OtherRelations &...other_relations)
    : UpdateRelation<ExecutionPolicy, FirstRelation>(first_relation),
      other_interactions_(other_relations...) {}
//=================================================================================================//
template <class ExecutionPolicy, class FirstRelation, class... OtherRelations>
void UpdateRelation<ExecutionPolicy, FirstRelation, OtherRelations...>::exec(Real dt)
{
    UpdateRelation<ExecutionPolicy, FirstRelation>::exec(dt);
    other_interactions_.exec(dt);
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_BODY_RELATION_HPP
