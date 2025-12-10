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
      src_pos_(encloser.inner_relation_.dvSourcePosition()->DelegatedData(ex_policy)),
      neighbor_size_(encloser.inner_relation_.dvNeighborSize()->DelegatedData(ex_policy)),
      masked_src_(ex_policy, encloser.inner_relation_.getDynamicsIdentifier()),
      is_one_sided_(ex_policy, encloser.inner_relation_.getNeighborhood()),
      masked_criterion_(
          ex_policy, encloser.inner_relation_.getDynamicsIdentifier(),
          ex_policy, encloser.inner_relation_.getNeighborhood()),
      neighbor_search_(encloser.cell_linked_list_.createNeighborSearch(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    InteractKernel::clearNeighborSize(UnsignedInt src_index)
{
    this->neighbor_index_[src_index] = 0;
    neighbor_size_[src_index] = 0;
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    InteractKernel::incrementNeighborSize(UnsignedInt src_index)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    if (masked_src_(src_index))
    {
        neighbor_search_.forEachSearch(
            src_pos_[src_index],
            [&](size_t tar_index)
            {
                if (is_one_sided_(src_index, tar_index) && masked_criterion_(tar_index, src_index))
                {
                    ++neighbor_count;
                    AtomicRef<UnsignedInt> atomic_tar_size(this->neighbor_index_[tar_index]);
                    ++atomic_tar_size;
                }
            });
    }
    AtomicRef<UnsignedInt> atomic_src_size(this->neighbor_index_[src_index]);
    atomic_src_size.fetch_add(neighbor_count);
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    InteractKernel::updateNeighborList(UnsignedInt src_index)
{
    if (masked_src_(src_index))
    {
        neighbor_search_.forEachSearch(
            src_pos_[src_index],
            [&](size_t tar_index)
            {
                if (is_one_sided_(src_index, tar_index) && masked_criterion_(tar_index, src_index))
                {
                    AtomicRef<UnsignedInt> atomic_src_size(this->neighbor_size_[src_index]);
                    this->neighbor_index_[this->particle_offset_[src_index] + atomic_src_size++] = tar_index;
                    AtomicRef<UnsignedInt> atomic_tar_size(this->neighbor_size_[tar_index]);
                    this->neighbor_index_[this->particle_offset_[tar_index] + atomic_tar_size++] = src_index;
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
                 { computing_kernel->clearNeighborSize(i); });

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementNeighborSize(i); });

    this->logger_->debug("UpdateCellLinkedList: incrementNeighborSize done at {} for Relation {}.",
                         this->sph_body_->getName(), type_name<Inner<Parameters...>>());

    auto *dv_neighbor_index = this->inner_relation_.dvNeighborIndex();
    auto *dv_particle_offset = this->inner_relation_.dvParticleOffset();
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
            old_size, dv_neighbor_index->getDataSize(), this->sph_body_->getName());
    }

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborList(i); });

    this->logger_->debug("UpdateCellLinkedList: updateNeighborList done at {} for Relation {}.",
                         this->sph_body_->getName(), type_name<Inner<Parameters...>>());
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
      src_pos_(encloser.contact_relation_.dvSourcePosition()->DelegatedData(ex_policy)),
      masked_src_(ex_policy, encloser.contact_relation_.getSourceIdentifier()),
      masked_criterion_(
          ex_policy, encloser.contact_relation_.getContactIdentifier(contact_index),
          ex_policy, encloser.contact_relation_.getNeighborhood(contact_index)),
      neighbor_search_(
          encloser.contact_cell_linked_list_[contact_index]->createNeighborSearch(ex_policy)),
      search_box_(ex_policy, encloser.contact_relation_.getNeighborhood(contact_index)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::incrementNeighborSize(UnsignedInt src_index)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    if (masked_src_(src_index))
    {
        neighbor_search_.forEachSearch(
            src_pos_[src_index],
            [&](size_t tar_index)
            {
                if (masked_criterion_(tar_index, src_index))
                    neighbor_count++;
            },
            search_box_(src_index));
    }
    this->neighbor_index_[src_index] = neighbor_count;
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::updateNeighborList(UnsignedInt src_index)
{
    UnsignedInt neighbor_count = 0;
    if (masked_src_(src_index))
    {
        neighbor_search_.forEachSearch(
            src_pos_[src_index],
            [&](size_t tar_index)
            {
                if (masked_criterion_(tar_index, src_index))
                {
                    this->neighbor_index_[this->particle_offset_[src_index] + neighbor_count] = tar_index;
                    neighbor_count++;
                }
            },
            search_box_(src_index));
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
                             this->sph_body_->getName(), type_name<Contact<Parameters...>>(),
                             contact_relation_.getContactIdentifier(k).getName());

        auto *dv_neighbor_index = this->contact_relation_.dvNeighborIndex(k);
        auto *dv_particle_offset = this->contact_relation_.dvParticleOffset(k);
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
                old_size, dv_neighbor_index->getDataSize(), this->sph_body_->getName());
        }

        particle_for(ex_policy_,
                     IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { computing_kernel->updateNeighborList(i); });

        this->logger_->debug(
            "UpdateRelation: updateNeighborList done at {} for Relation {} to {}.",
            this->sph_body_->getName(), type_name<Contact<Parameters...>>(),
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
