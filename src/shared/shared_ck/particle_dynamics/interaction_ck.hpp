#ifndef INTERACTION_CK_HPP
#define INTERACTION_CK_HPP

#include "interaction_ck.h"

#include "base_material.h"

namespace SPH
{
//=================================================================================================//
template <typename... Parameters>
Interaction<Inner<Parameters...>>::
    Interaction(InnerRelationType &inner_relation)
    : BaseLocalDynamicsType(inner_relation.getDynamicsIdentifier()),
      inner_relation_(&inner_relation),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Inner<Parameters...>>::
    registerComputingKernel(Implementation<Base> *implementation)
{
    inner_relation_->registerComputingKernel(implementation);
}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Inner<Parameters...>>::resetComputingKernelUpdated()
{
    inner_relation_->resetComputingKernelUpdated();
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interaction<Inner<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : NeighborList(ex_policy, *encloser.inner_relation_),
      NeighborKernel(ex_policy, encloser.inner_relation_->getNeighborhood()) {}
//=================================================================================================//
template <typename... Parameters>
Interaction<Contact<Parameters...>>::Interaction(Contact<Parameters...> &contact_relation)
    : BaseLocalDynamicsType(contact_relation.getSourceIdentifier()),
      contact_relation_(&contact_relation),
      contact_body_(&contact_relation.getContactBody()),
      contact_particles_(&contact_relation.getContactParticles()),
      contact_adaptation_(&contact_relation.getContactAdaptation()),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_contact_Vol_(
          contact_particles_->template getVariableByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Contact<Parameters...>>::registerComputingKernel(Implementation<Base> *implementation)
{
    contact_relation_->registerComputingKernel(implementation);
}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Contact<Parameters...>>::resetComputingKernelUpdated()
{
    contact_relation_->resetComputingKernelUpdated();
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interaction<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : NeighborList(ex_policy, *encloser.contact_relation_),
      NeighborKernel(ex_policy, encloser.contact_relation_->getNeighborhood()) {}
//=================================================================================================//
template <class WallContactRelationType>
Interaction<Wall>::Interaction(WallContactRelationType &wall_contact_relation)
{
    BaseParticles &contact_particles = wall_contact_relation.getContactParticles();
    SPHBody &contact_body = wall_contact_relation.getContactBody();
    Solid &solid_material = DynamicCast<Solid>(this, contact_body.getMatterMaterial());
    dv_wall_vel_ave_ = solid_material.AverageVelocityVariable(&contact_particles);
    dv_wall_acc_ave_ = solid_material.AverageAccelerationVariable(&contact_particles);
    dv_wall_n_ = contact_particles.template getVariableByName<Vecd>("NormalDirection");
}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_CK_HPP