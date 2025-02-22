#ifndef RELATION_CK_HPP
#define RELATION_CK_HPP

#include "relation_ck.h"

namespace SPH
{
//=================================================================================================//
template <class DataType>
DiscreteVariable<DataType> *Relation<Base>::
    addRelationVariable(const std::string &name, size_t data_size)
{
    return relation_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
}
//=================================================================================================//
template <class DynamicsIdentifier, class TargetIdentifier>
Relation<Contact<DynamicsIdentifier, TargetIdentifier>>::
    Relation(DynamicsIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers)
    : Relation<Base>(source_identifier.getSPHBody()),
      source_identifier_(source_identifier), contact_identifiers_(contact_identifiers)
{
    for (size_t k = 0; k != contact_identifiers.size(); ++k)
    {
        RealBody *contact_body = DynamicCast<RealBody>(this, &contact_identifiers[k]->getSPHBody());

        contact_bodies_.push_back(contact_body);
        const std::string name = contact_body->getName();
        contact_particles_.push_back(&contact_body->getBaseParticles());
        contact_adaptations_.push_back(&contact_body->getSPHAdaptation());

        CellLinkedList *target_cell_linked_list =
            DynamicCast<CellLinkedList>(this, &contact_body->getCellLinkedList());
        target_cell_linked_lists_.push_back(target_cell_linked_list);

        dv_contact_neighbor_index_.push_back(addRelationVariable<UnsignedInt>(
            "Contact" + name + "NeighborIndex", offset_list_size_));
        dv_contact_particle_offset_.push_back(addRelationVariable<UnsignedInt>(
            "Contact" + name + "ParticleOffset", offset_list_size_));
        all_contact_computing_kernels_.resize(contact_bodies_.size());
    }
}
//=================================================================================================//
template <class DynamicsIdentifier, class TargetIdentifier>
void Relation<Contact<DynamicsIdentifier, TargetIdentifier>>::registerComputingKernel(
    execution::Implementation<Base> *implementation, UnsignedInt contact_index)
{
    all_contact_computing_kernels_[contact_index].push_back(implementation);
}
//=================================================================================================//
template <class DynamicsIdentifier, class TargetIdentifier>
void Relation<Contact<DynamicsIdentifier, TargetIdentifier>>::
    resetComputingKernelUpdated(UnsignedInt contact_index)
{
    for (size_t k = 0; k != all_contact_computing_kernels_[contact_index].size(); ++k)
    {
        all_contact_computing_kernels_[contact_index][k]->resetUpdated();
    }
}
//=================================================================================================//
} // namespace SPH
#endif // RELATION_CK_HPP
