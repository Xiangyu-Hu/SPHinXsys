#ifndef RELATION_CK_HPP
#define RELATION_CK_HPP

#include "relation_ck.h"

namespace SPH
{
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier>
Relation<SourceIdentifier, TargetIdentifier>::Relation(
    SourceIdentifier &src_identifier, TargetIdentifier &tgt_identifier, ConfigType config_type)
    : RelationBase(src_identifier.Name(), tgt_identifier.Name()),
      sph_body_(&src_identifier.getSPHBody()), particles_(&sph_body_->getBaseParticles()),
      dv_source_pos_(this->assignConfigPosition(*particles_, config_type)),
      dv_neighbor_size_(particles_->registerDiscreteVariable<UnsignedInt>(
          "NeighborSize", particles_->ParticlesBound())),
      offset_list_size_(particles_->ParticlesBound() + 1)
{
    SPHBody &tgt_body = tgt_identifier.getSPHBody();
    BaseParticles &tgt_particles = tgt_body.getBaseParticles();
    dv_target_pos_ = assignConfigPosition(tgt_particles, config_type);
    dv_target_neighbor_index_ = addRelationVariable<UnsignedInt>(
        name_ + "NeighborIndex", offset_list_size_);
    dv_target_particle_offset_ = addRelationVariable<UnsignedInt>(
        name_ + "ParticleOffset", offset_list_size_);
    neighborhood_ = neighborhood_ptr_.template createPtr<NeighborhoodType>(
        src_identifier, tgt_identifier, dv_source_pos_, dv_target_pos_);
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier>
DiscreteVariable<Vecd> *Relation<SourceIdentifier, TargetIdentifier>::
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
template <class SourceIdentifier, class TargetIdentifier>
template <class DataType>
DiscreteVariable<DataType> *Relation<SourceIdentifier, TargetIdentifier>::
    addRelationVariable(const std::string &name, size_t data_size)
{
    return relation_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier>
void Relation<SourceIdentifier, TargetIdentifier>::registerComputingKernel(
    execution::Implementation<Base> *implementation)
{
    registered_computing_kernels_.push_back(implementation);
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier>
void Relation<SourceIdentifier, TargetIdentifier>::resetComputingKernelUpdated()
{
    for (size_t k = 0; k != registered_computing_kernels_.size(); ++k)
    {
        registered_computing_kernels_[k]->resetUpdated();
    }
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier>
template <class ExecutionPolicy, class EncloserType>
Relation<SourceIdentifier, TargetIdentifier>::NeighborList::NeighborList(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : neighbor_index_(encloser.dv_target_neighbor_index_->DelegatedData(ex_policy)),
      particle_offset_(encloser.dv_target_particle_offset_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DynamicsIdentifier>
template <typename... Args>
Inner<Relation<DynamicsIdentifier>>::
    Inner(DynamicsIdentifier &identifier, Args &&...args)
    : Relation<DynamicsIdentifier, DynamicsIdentifier>(
          identifier, identifier, std::forward<Args>(args)...),
      identifier_(&identifier) {}
//=================================================================================================//
template <typename SourceIdentifier, class TargetIdentifier>
Contact<Relation<SourceIdentifier, TargetIdentifier>>::Contact(
    SourceIdentifier &src_identifier, TargetIdentifier &tgt_identifier, ConfigType config_type)
    : Relation<SourceIdentifier, TargetIdentifier>(src_identifier, tgt_identifier, config_type),
      src_identifier_(&src_identifier), tgt_identifier_(&tgt_identifier),
      contact_body_(&tgt_identifier.getSPHBody()),
      contact_particles_(&contact_body_->getBaseParticles()),
      contact_adaptation_(&contact_body_->getSPHAdaptation()) {}
//=================================================================================================//
} // namespace SPH
#endif // RELATION_CK_HPP
