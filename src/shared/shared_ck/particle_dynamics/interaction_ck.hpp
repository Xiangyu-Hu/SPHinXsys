#ifndef INTERACTION_CK_HPP
#define INTERACTION_CK_HPP

#include "interaction_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename... Parameters>
Interaction<Inner<Parameters...>>::
    Interaction(Relation<Inner<Parameters...>> &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      inner_relation_(inner_relation),
      real_body_(&inner_relation.getRealBody()),
      sph_adaptation_(&sph_body_.getSPHAdaptation()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_neighbor_index_(inner_relation.getNeighborIndex()),
      dv_particle_offset_(inner_relation.getParticleOffset()) {}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Inner<Parameters...>>::
    registerComputingKernel(Implementation<Base> *implementation)
{
    inner_relation_.registerComputingKernel(implementation);
}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Inner<Parameters...>>::resetComputingKernelUpdated()
{
    inner_relation_.resetComputingKernelUpdated();
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
Interaction<Inner<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   Interaction<Inner<Parameters...>> &encloser)
    : NeighborList(ex_policy, encloser.dv_neighbor_index_, encloser.dv_particle_offset_),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_, encloser.dv_pos_) {}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Parameters>
Interaction<Contact<SourceIdentifier, TargetIdentifier, Parameters...>>::
    Interaction(Relation<Contact<SourceIdentifier, TargetIdentifier, Parameters...>> &contact_relation)
    : BaseLocalDynamics<SourceIdentifier>(contact_relation.getSourceIdentifier()),
      contact_relation_(contact_relation),
      sph_adaptation_(&this->sph_body_.getSPHAdaptation()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      contact_bodies_(contact_relation.getContactBodies()),
      contact_particles_(contact_relation.getContactParticles()),
      contact_adaptations_(contact_relation.getContactAdaptations()),
      dv_contact_neighbor_index_(contact_relation.getTargetNeighborIndex()),
      dv_contact_particle_offset_(contact_relation.getTargetParticleOffset())
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_pos_.push_back(contact_particles_[k]->template getVariableByName<Vecd>("Position"));
    }
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Parameters>
void Interaction<Contact<SourceIdentifier, TargetIdentifier, Parameters...>>::
    registerComputingKernel(Implementation<Base> *implementation, UnsignedInt contact_index)
{
    contact_relation_.registerComputingKernel(implementation, contact_index);
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Parameters>
void Interaction<Contact<SourceIdentifier, TargetIdentifier, Parameters...>>::
    resetComputingKernelUpdated(UnsignedInt contact_index)
{
    contact_relation_.resetComputingKernelUpdated(contact_index);
}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interaction<Contact<SourceIdentifier, TargetIdentifier, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : NeighborList(ex_policy, encloser.dv_contact_neighbor_index_[contact_index],
                   encloser.dv_contact_particle_offset_[contact_index]),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_,
                              encloser.contact_adaptations_[contact_index],
                              encloser.dv_pos_, encloser.contact_pos_[contact_index]) {}
//=================================================================================================//
template <typename... Parameters>
Interaction<Contact<Wall, Parameters...>>::
    Interaction(Relation<Contact<Parameters...>> &wall_contact_relation)
    : Interaction<Contact<Parameters...>>(wall_contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Solid &solid_material = DynamicCast<Solid>(this, this->contact_particles_[k]->getBaseMaterial());
        dv_wall_vel_ave_.push_back(solid_material.AverageVelocityVariable(this->contact_particles_[k]));
        dv_wall_acc_ave_.push_back(solid_material.AverageAccelerationVariable(this->contact_particles_[k]));
        dv_wall_n_.push_back(this->contact_particles_[k]->template getVariableByName<Vecd>("NormalDirection"));
        dv_wall_Vol_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_CK_HPP
