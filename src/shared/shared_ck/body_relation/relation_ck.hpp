#ifndef RELATION_CK_HPP
#define RELATION_CK_HPP

#include "relation_ck.h"

namespace SPH
{
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier>
Relation<Base>::Relation(SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers)
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
template <class DynamicsIdentifier, class TargetIdentifier>
Relation<Contact<DynamicsIdentifier, TargetIdentifier>>::
    Relation(DynamicsIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers)
    : Relation<Base>(source_identifier, contact_identifiers),
      source_identifier_(source_identifier), contact_identifiers_(contact_identifiers),
      dv_contact_neighbor_index_(dv_target_neighbor_index_),
      dv_contact_particle_offset_(dv_target_particle_offset_)
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
template <class DynamicsIdentifier, class TargetIdentifier>
void Relation<Contact<DynamicsIdentifier, TargetIdentifier>>::registerComputingKernel(
    execution::Implementation<Base> *implementation, UnsignedInt contact_index)
{
    registered_computing_kernels_[contact_index].push_back(implementation);
}
//=================================================================================================//
template <class DynamicsIdentifier, class TargetIdentifier>
void Relation<Contact<DynamicsIdentifier, TargetIdentifier>>::
    resetComputingKernelUpdated(UnsignedInt contact_index)
{
    for (size_t k = 0; k != registered_computing_kernels_[contact_index].size(); ++k)
    {
        registered_computing_kernels_[contact_index][k]->resetUpdated();
    }
}
//=================================================================================================//
} // namespace SPH
#endif // RELATION_CK_HPP
