#ifndef INTERACTION_CK_HPP
#define INTERACTION_CK_HPP

#include "interaction_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename... Parameters>
Interaction<Inner<Parameters...>>::
    Interaction(InnerRelationType &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      inner_relation_(inner_relation),
      real_body_(&inner_relation.getRealBody()) {}
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
    InteractKernel(const ExecutionPolicy &ex_policy, Interaction<Inner<Parameters...>> &encloser)
    : NeighborList(ex_policy, encloser.inner_relation_),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_,
                              encloser.sph_adaptation_,
                              encloser.inner_relation_.getSourcePosition(),
                              encloser.inner_relation_.getTargetPosition()) {}
//=================================================================================================//
template <class SourceIdentifier, class TargetIdentifier, typename... Parameters>
Interaction<Contact<SourceIdentifier, TargetIdentifier, Parameters...>>::
    Interaction(ContactRelationType &contact_relation)
    : BaseLocalDynamics<SourceIdentifier>(contact_relation.getSourceIdentifier()),
      contact_relation_(contact_relation),
      sph_adaptation_(&this->sph_body_.getSPHAdaptation()),
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
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : NeighborList(ex_policy, encloser.contact_relation_, contact_index),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_,
                              encloser.contact_adaptations_[contact_index],
                              encloser.contact_relation_.getSourcePosition(),
                              encloser.contact_relation_.getTargetPosition(contact_index)) {}
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
