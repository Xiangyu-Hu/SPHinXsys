#ifndef INTERACTION_CK_HPP
#define INTERACTION_CK_HPP

#include "interaction_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename... Parameters>
Interaction<Inner<Parameters...>>::
    Interaction(InnerRelationType &inner_relation)
    : BaseLocalDynamicsType(inner_relation.getDynamicsIdentifier()),
      inner_relation_(inner_relation),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
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
template <class ExecutionPolicy, class EncloserType>
Interaction<Inner<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : NeighborList(ex_policy, encloser.inner_relation_),
      NeighborKernel(ex_policy, encloser.inner_relation_.getNeighborhood()) {}
//=================================================================================================//
template <typename... Parameters>
Interaction<Contact<Parameters...>>::
    Interaction(ContactRelationType &contact_relation)
    : BaseLocalDynamicsType(contact_relation.getSourceIdentifier()),
      contact_relation_(contact_relation),
      contact_bodies_(contact_relation.getContactBodies()),
      contact_particles_(contact_relation.getContactParticles()),
      contact_adaptations_(contact_relation.getContactAdaptations()),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure"))
{
    for (auto &particles : contact_particles_)
    {
        dv_contact_Vol_.push_back(
            particles->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Contact<Parameters...>>::
    registerComputingKernel(Implementation<Base> *implementation, UnsignedInt contact_index)
{
    contact_relation_.registerComputingKernel(implementation, contact_index);
}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Contact<Parameters...>>::
    resetComputingKernelUpdated(UnsignedInt contact_index)
{
    contact_relation_.resetComputingKernelUpdated(contact_index);
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interaction<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : NeighborList(ex_policy, encloser.contact_relation_, contact_index),
      NeighborKernel(ex_policy, encloser.contact_relation_.getNeighborhood(contact_index)) {}
//=================================================================================================//
template <class WallContactRelationType>
Interaction<Wall>::Interaction(WallContactRelationType &wall_contact_relation)
{
    StdVec<BaseParticles *> contact_particles = wall_contact_relation.getContactParticles();
    for (size_t k = 0; k != contact_particles.size(); ++k)
    {
        Solid &solid_material = DynamicCast<Solid>(this, contact_particles[k]->getBaseMaterial());
        dv_wall_vel_ave_.push_back(solid_material.AverageVelocityVariable(contact_particles[k]));
        dv_wall_acc_ave_.push_back(solid_material.AverageAccelerationVariable(contact_particles[k]));
        dv_wall_n_.push_back(contact_particles[k]->template getVariableByName<Vecd>("NormalDirection"));
    }
}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_CK_HPP
