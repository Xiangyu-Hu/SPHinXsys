#ifndef RELATION_CK_HPP
#define RELATION_CK_HPP

#include "relation_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename... AdaptationParameters>
template <class SourceIdentifier, class TargetIdentifier>
Relation<NeighborMethod<AdaptationParameters...>>::Relation(
    SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers, ConfigType config_type)
    : sph_body_(&source_identifier.getSPHBody()),
      particles_(&sph_body_->getBaseParticles()),
      dv_source_pos_(this->assignConfigPosition(*particles_, config_type)),
      dv_neighbor_size_(addRelationVariable<UnsignedInt>(
          source_identifier.getName() + "NeighborSize", particles_->ParticlesBound())),
      offset_list_size_(particles_->ParticlesBound() + 1)
{
    std::string source_name = source_identifier.getName();
    for (size_t k = 0; k != contact_identifiers.size(); ++k)
    {
        SPHBody &contact_body = contact_identifiers[k]->getSPHBody();
        std::string target_name = contact_identifiers[k]->getName();
        std::string name = source_name == target_name
                               ? source_name + "Inner"
                               : source_name + "To" + target_name;
        BaseParticles &contact_particles = contact_body.getBaseParticles();
        dv_target_pos_.push_back(assignConfigPosition(contact_particles, config_type));
        dv_target_neighbor_index_.push_back(addRelationVariable<UnsignedInt>(
            name + "NeighborIndex", offset_list_size_));
        dv_target_particle_offset_.push_back(addRelationVariable<UnsignedInt>(
            name + "ParticleOffset", offset_list_size_));
        neighborhoods_.push_back(
            neighborhood_ptrs_.template createPtr<Neighbor<NeighborMethodType>>(
                source_identifier, *contact_identifiers[k], dv_source_pos_, dv_target_pos_.back()));
    }
    registered_computing_kernels_.resize(contact_identifiers.size());
}
//=================================================================================================//
template <typename... AdaptationParameters>
DiscreteVariable<Vecd> *Relation<NeighborMethod<AdaptationParameters...>>::
    assignConfigPosition(BaseParticles &particles, ConfigType config_type)
{
    if (config_type == ConfigType::Eulerian)
    {
        return particles.getVariableByName<Vecd>("Position");
    }
    else
    {
        return particles.registerStateVariableFrom<Vecd>(
            "InitialPosition", "Position");
    }
}
//=================================================================================================//
template <typename... AdaptationParameters>
template <class DataType>
DiscreteVariable<DataType> *Relation<NeighborMethod<AdaptationParameters...>>::
    addRelationVariable(const std::string &name, size_t data_size)
{
    return relation_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
}
//=================================================================================================//
template <typename... AdaptationParameters>
void Relation<NeighborMethod<AdaptationParameters...>>::registerComputingKernel(
    execution::Implementation<Base> *implementation, UnsignedInt target_index)
{
    registered_computing_kernels_[target_index].push_back(implementation);
}
//=================================================================================================//
template <typename... AdaptationParameters>
void Relation<NeighborMethod<AdaptationParameters...>>::
    resetComputingKernelUpdated(UnsignedInt target_index)
{
    auto &computing_kernels = registered_computing_kernels_[target_index];
    for (size_t k = 0; k != computing_kernels.size(); ++k)
    {
        computing_kernels[k]->resetUpdated();
    }
}
//=================================================================================================//
template <typename... AdaptationParameters>
template <class ExecutionPolicy, class EncloserType>
Relation<NeighborMethod<AdaptationParameters...>>::NeighborList::NeighborList(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt target_index)
    : neighbor_index_(encloser.dv_target_neighbor_index_[target_index]->DelegatedData(ex_policy)),
      particle_offset_(encloser.dv_target_particle_offset_[target_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DynamicsIdentifier, typename... AdaptationParameters>
template <typename... Args>
Inner<DynamicsIdentifier, NeighborMethod<AdaptationParameters...>>::
    Inner(DynamicsIdentifier &identifier, Args &&...args)
    : Relation<NeighborMethod<AdaptationParameters...>>(
          identifier, StdVec<DynamicsIdentifier *>{&identifier}, std::forward<Args>(args)...),
      identifier_(&identifier) {}
//=================================================================================================//
template <typename SourceIdentifier, class TargetIdentifier, typename... AdaptationParameters>
Contact<SourceIdentifier, TargetIdentifier, NeighborMethod<AdaptationParameters...>>::Contact(
    SourceIdentifier &source_identifier, StdVec<TargetIdentifier *> contact_identifiers, ConfigType config_type)
    : Relation<NeighborMethod<AdaptationParameters...>>(source_identifier, contact_identifiers, config_type),
      source_identifier_(&source_identifier), contact_identifiers_(contact_identifiers)
{
    for (size_t k = 0; k != contact_identifiers.size(); ++k)
    {
        SPHBody *contact_body = &contact_identifiers[k]->getSPHBody();
        contact_bodies_.push_back(contact_body);
        contact_particles_.push_back(&contact_body->getBaseParticles());
        contact_adaptations_.push_back(&contact_body->getSPHAdaptation());
    }
}
//=================================================================================================//
} // namespace SPH
#endif // RELATION_CK_HPP
