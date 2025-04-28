#ifndef INTERACTION_CK_HPP
#define INTERACTION_CK_HPP

#include "interaction_ck.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsIdentifier, typename... Parameters>
Interaction<Inner<DynamicsIdentifier, Parameters...>>::
    Interaction(InnerRelationType &inner_relation)
    : BaseLocalDynamics<DynamicsIdentifier>(inner_relation.getDynamicsIdentifier()),
      inner_relation_(inner_relation) {}
//=================================================================================================//
template <class DynamicsIdentifier, typename... Parameters>
void Interaction<Inner<DynamicsIdentifier, Parameters...>>::
    registerComputingKernel(Implementation<Base> *implementation)
{
    inner_relation_.registerComputingKernel(implementation);
}
//=================================================================================================//
template <class DynamicsIdentifier, typename... Parameters>
void Interaction<Inner<DynamicsIdentifier, Parameters...>>::resetComputingKernelUpdated()
{
    inner_relation_.resetComputingKernelUpdated();
}
//=================================================================================================//
template <class DynamicsIdentifier, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interaction<Inner<DynamicsIdentifier, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Parameters &&...parameters)
    : NeighborList(ex_policy, encloser.inner_relation_),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_,
                              encloser.sph_adaptation_,
                              encloser.inner_relation_.getSourcePosition(),
                              encloser.inner_relation_.getTargetPosition(),
                              std::forward<Parameters>(parameters)...) {}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Parameters>
Interaction<Contact<SourceIdentifier, TargetIdentifier, Parameters...>>::
    Interaction(ContactRelationType &contact_relation)
    : BaseLocalDynamics<SourceIdentifier>(contact_relation.getSourceIdentifier()),
      contact_relation_(contact_relation),
      contact_bodies_(contact_relation.getContactBodies()),
      contact_particles_(contact_relation.getContactParticles()),
      contact_adaptations_(contact_relation.getContactAdaptations()) {}
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
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                   UnsignedInt contact_index, Parameters &&...parameters)
    : NeighborList(ex_policy, encloser.contact_relation_, contact_index),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_,
                              encloser.contact_adaptations_[contact_index],
                              encloser.contact_relation_.getSourcePosition(),
                              encloser.contact_relation_.getTargetPosition(contact_index),
                              std::forward<Parameters>(parameters)...) {}
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
