#ifndef UPDATE_BODY_RELATION_HPP
#define UPDATE_BODY_RELATION_HPP

#include "update_body_relation.h"

#include "cell_linked_list.hpp"
#include "particle_iterators_ck.h"
#include "sphinxsys_atom_ref.h"

namespace SPH
{
//=================================================================================================//
template <typename TargetCriterion>
TargetParticleMask<TargetCriterion, BodyPartByParticle>::TargetParticleMask(
    BodyPartByParticle &body_by_particle)
    : particle_group_manager_(body_by_particle.getParticleGroupManager()),
      body_part_bit_(particle_group_manager_.getGroupMask(body_by_particle.PartIDName())) {}
//=================================================================================================//
template <typename TargetCriterion>
template <class ExecutionPolicy, typename EnclosureType, typename... Args>
TargetParticleMask<TargetCriterion, BodyPartByParticle>::ComputingKernel::
    ComputingKernel(ExecutionPolicy &ex_policy, EnclosureType &encloser, Args &&...args)
    : TargetCriterion(std::forward<Args>(args)...),
      body_part_mask_(ex_policy, encloser.particle_group_manager_, encloser.body_part_bit_) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
UpdateRelation<ExecutionPolicy, Inner<Parameters...>>::
    UpdateRelation(Inner<Parameters...> &inner_relation)
    : BaseLocalDynamicsType(inner_relation.getDynamicsIdentifier()),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{}),
      inner_relation_(inner_relation),
      target_particle_mask_method_(inner_relation_.getDynamicsIdentifier()),
      cell_linked_list_(DynamicCast<CellLinkedList<CellLinkedListIdentifier>>(
          this, inner_relation_.getDynamicsIdentifier().getCellLinkedList())),
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
          ex_policy, encloser.target_particle_mask_method_,
          ex_policy, encloser.inner_relation_.getNeighborhood()),
      neighbor_search_(ex_policy, encloser.cell_linked_list_),
      src_cut_off_(ex_policy, encloser.inner_relation_.getNeighborhood()) {}
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
        neighbor_search_.forInnerSearch(
            src_pos_[src_index],
            [&](size_t tar_index)
            {
                if (is_one_sided_(src_index, tar_index) && masked_criterion_(tar_index, src_index))
                {
                    ++neighbor_count;
                    AtomicRef<UnsignedInt> atomic_tar_size(this->neighbor_index_[tar_index]);
                    ++atomic_tar_size;
                }
            },
            src_cut_off_(src_index));
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
        neighbor_search_.forInnerSearch(
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
            },
            src_cut_off_(src_index));
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
                         this->sph_body_->Name(), type_name<Inner<Parameters...>>());

    auto *dv_neighbor_index = this->inner_relation_.dvNeighborIndex();
    auto *dv_particle_offset = this->inner_relation_.dvParticleOffset();
    UnsignedInt *neighbor_index = dv_neighbor_index->DelegatedData(ex_policy_);
    UnsignedInt *particle_offset = dv_particle_offset->DelegatedData(ex_policy_);
    UnsignedInt current_offset_list_size = total_real_particles + 1;
    UnsignedInt current_neighbor_index_size =
        exclusive_scan(ex_policy_, neighbor_index, particle_offset, current_offset_list_size,
                       typename PlusUnsignedInt<ExecutionPolicy>::type());

    if (current_neighbor_index_size > dv_neighbor_index->getSize())
    {
        UnsignedInt old_size = dv_neighbor_index->getSize();
        dv_neighbor_index->reallocateData(ex_policy_, current_neighbor_index_size);
        this->inner_relation_.resetComputingKernelUpdated();
        kernel_implementation_.overwriteComputingKernel();

        this->logger_->info(
            "UpdateRelation: increase neighbor index size from {} to {} at {} .",
            old_size, dv_neighbor_index->getSize(), this->sph_body_->Name());
    }

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborList(i); });

    this->logger_->debug("UpdateCellLinkedList: updateNeighborList done at {} for Relation {}.",
                         this->sph_body_->Name(), type_name<Inner<Parameters...>>());
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    UpdateRelation(ContactRelationType &contact_relation)
    : BaseLocalDynamicsType(contact_relation.getSourceIdentifier()),
      BaseDynamics<void>(), ex_policy_(ExecutionPolicy{}),
      contact_relation_(contact_relation),
      contact_target_particle_mask_method_(contact_relation.getTargetIdentifier()),
      contact_cell_linked_list_(
          DynamicCast<CellLinkedList<CellLinkedListIdentifier>>(
              this, &contact_relation.getTargetIdentifier().getCellLinkedList())),
      contact_kernel_implementation_(
          contact_kernel_implementation_ptr_.template createPtr<KernelImplementation>(*this)) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
template <class EncloserType>
UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : NeighborList(ex_policy, encloser.contact_relation_),
      src_pos_(encloser.contact_relation_.dvSourcePosition()->DelegatedData(ex_policy)),
      masked_src_(ex_policy, encloser.contact_relation_.getSourceIdentifier()),
      masked_criterion_(
          ex_policy, encloser.contact_target_particle_mask_method_,
          ex_policy, encloser.contact_relation_.getNeighborhood()),
      neighbor_search_(ex_policy, *encloser.contact_cell_linked_list_),
      src_cut_off_(ex_policy, encloser.contact_relation_.getNeighborhood()) {}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::
    InteractKernel::incrementNeighborSize(UnsignedInt src_index)
{
    // Here, neighbor_index_ takes role of temporary storage for neighbor size list.
    UnsignedInt neighbor_count = 0;
    if (masked_src_(src_index))
    {
        neighbor_search_.forContactSearch(
            src_pos_[src_index],
            [&](size_t tar_index)
            {
                if (masked_criterion_(tar_index, src_index))
                    neighbor_count++;
            },
            src_cut_off_(src_index));
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
        neighbor_search_.forContactSearch(
            src_pos_[src_index],
            [&](size_t tar_index)
            {
                if (masked_criterion_(tar_index, src_index))
                {
                    this->neighbor_index_[this->particle_offset_[src_index] + neighbor_count] = tar_index;
                    neighbor_count++;
                }
            },
            src_cut_off_(src_index));
    }
}
//=================================================================================================//
template <class ExecutionPolicy, typename... Parameters>
void UpdateRelation<ExecutionPolicy, Contact<Parameters...>>::exec(Real dt)
{
    UnsignedInt total_real_particles = this->particles_->TotalRealParticles();

    InteractKernel *computing_kernel = contact_kernel_implementation_->getComputingKernel();
    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementNeighborSize(i); });

    this->logger_->debug("UpdateRelation: incrementNeighborSize done at {} for Relation {} to {}.",
                         this->sph_body_->Name(), type_name<Contact<Parameters...>>(),
                         contact_relation_.getTargetIdentifier().Name());

    auto *dv_neighbor_index = this->contact_relation_.dvNeighborIndex();
    auto *dv_particle_offset = this->contact_relation_.dvParticleOffset();
    UnsignedInt *neighbor_index = dv_neighbor_index->DelegatedData(ex_policy_);
    UnsignedInt *particle_offset = dv_particle_offset->DelegatedData(ex_policy_);
    UnsignedInt current_offset_list_size = total_real_particles + 1;
    UnsignedInt current_neighbor_index_size =
        exclusive_scan(ex_policy_, neighbor_index, particle_offset, current_offset_list_size,
                       typename PlusUnsignedInt<ExecutionPolicy>::type());

    if (current_neighbor_index_size > dv_neighbor_index->getSize())
    {
        UnsignedInt old_size = dv_neighbor_index->getSize();
        dv_neighbor_index->reallocateData(ex_policy_, current_neighbor_index_size);
        this->contact_relation_.resetComputingKernelUpdated();
        contact_kernel_implementation_->overwriteComputingKernel();

        this->logger_->info(
            "UpdateRelation: increase neighbor index size from {} to {} at .",
            old_size, dv_neighbor_index->getSize(), this->sph_body_->Name());
    }

    particle_for(ex_policy_,
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateNeighborList(i); });

    this->logger_->debug(
        "UpdateRelation: updateNeighborList done at {} for Relation {} to {}.",
        this->sph_body_->Name(), type_name<Contact<Parameters...>>(),
        contact_relation_.getTargetIdentifier().Name());
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
