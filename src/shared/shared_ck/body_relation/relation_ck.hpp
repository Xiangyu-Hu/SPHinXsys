#ifndef RELATION_CK_HPP
#define RELATION_CK_HPP

#include "relation_ck.h"

namespace SPH
{
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier>
Relation<Base>::Relation(
    SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers)
    : sph_body_(source_identifier.getSPHBody()),
      particles_(sph_body_.getBaseParticles()),
      offset_list_size_(particles_.ParticlesBound() + 1)
{
    for (size_t k = 0; k != contact_identifiers.size(); ++k)
    {
        SPHBody &contact_body = contact_identifiers[k]->getSPHBody();
        const std::string name = contact_body.getName();
        dv_target_neighbor_index_.push_back(addRelationVariable<UnsignedInt>(
            name + "NeighborIndex", offset_list_size_));
        dv_target_particle_offset_.push_back(addRelationVariable<UnsignedInt>(
            name + "ParticleOffset", offset_list_size_));
    }
    registered_computing_kernels_.resize(contact_identifiers.size());
}
//=================================================================================================//
template <class DataType>
DiscreteVariable<DataType> *Relation<Base>::
    addRelationVariable(const std::string &name, size_t data_size)
{
    return relation_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
Relation<Base>::NeighborList::NeighborList(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt target_index)
    : neighbor_index_(encloser.dv_target_neighbor_index_[target_index]->DelegatedData(ex_policy)),
      particle_offset_(encloser.dv_target_particle_offset_[0]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier, class TargetIdentifier>
Relation<Contact<DynamicsIdentifier, TargetIdentifier>>::Relation(
    DynamicsIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers)
    : Relation<Base>(source_identifier, contact_identifiers),
      source_identifier_(source_identifier), contact_identifiers_(contact_identifiers)
{
    for (size_t k = 0; k != contact_identifiers.size(); ++k)
    {
        RealBody *contact_body = DynamicCast<RealBody>(this, &contact_identifiers[k]->getSPHBody());
        contact_bodies_.push_back(contact_body);
        contact_particles_.push_back(&contact_body->getBaseParticles());
        contact_adaptations_.push_back(&contact_body->getSPHAdaptation());
    }
}
//=================================================================================================//
} // namespace SPH
#endif // RELATION_CK_HPP
