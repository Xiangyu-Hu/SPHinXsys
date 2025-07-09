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
      inner_relation_(inner_relation) {}
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
      Neighborhood(ex_policy, encloser.inner_relation_.getNeighborMethod()) {}
//=================================================================================================//
template <typename... Parameters>
Interaction<Contact<Parameters...>>::
    Interaction(ContactRelationType &contact_relation)
    : BaseLocalDynamicsType(contact_relation.getSourceIdentifier()),
      contact_relation_(contact_relation),
      contact_bodies_(contact_relation.getContactBodies()),
      contact_particles_(contact_relation.getContactParticles()),
      contact_adaptations_(contact_relation.getContactAdaptations()) {}
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
      Neighborhood(ex_policy, encloser.contact_relation_.getNeighborMethod(contact_index)) {}
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
        dv_wall_Vol_.push_back(contact_particles[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <class SoilContactRelationType>
Interaction<Soil>::Interaction(SoilContactRelationType &soil_contact_relation)
{
    StdVec<BaseParticles *> contact_particles = soil_contact_relation.getContactParticles();
    for (size_t k = 0; k != contact_particles.size(); ++k)
    {
        dv_soil_n_.push_back(contact_particles[k]->template getVariableByName<Vecd>("NormalDirection"));
        dv_soil_vel_.push_back(contact_particles[k]->template getVariableByName<Vecd>("Velocity"));
        dv_soil_Vol_.push_back(contact_particles[k]->template getVariableByName<Real>("VolumetricMeasure"));
        dv_soil_p_.push_back(contact_particles[k]->template getVariableByName<Real>("Pressure"));
    }
}
//=================================================================================================//
template <class FluidContactRelationType>
Interaction<Fluid>::Interaction(FluidContactRelationType &fluid_contact_relation)
{
    StdVec<BaseParticles *> contact_particles = fluid_contact_relation.getContactParticles();
    for (size_t k = 0; k != contact_particles.size(); ++k)
    {
        dv_fluid_vel_.push_back(contact_particles[k]->template registerStateVariableOnly<Vecd>("Velocity"));
        dv_fluid_Vol_.push_back(contact_particles[k]->template registerStateVariableOnly<Real>("VolumetricMeasure"));
        dv_fluid_p_.push_back(contact_particles[k]->template registerStateVariableOnly<Real>("Pressure"));
    }
}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_CK_HPP
