#ifndef INTERACTION_CK_HPP
#define INTERACTION_CK_HPP

#include "interaction_ck.h"

#include "base_material.h"

#include <cstdlib>
#include <iostream>

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
Interaction<Contact<Parameters...>>::Interaction(ContactRelationType &contact_relation)
    : BaseLocalDynamicsType(contact_relation.getSourceIdentifier()),
      contact_relation_(&contact_relation),
      contact_bodies_(contact_relation.getContactBodies()),
      contact_particles_(contact_relation.getContactParticles()),
      contact_adaptations_(contact_relation.getContactAdaptations()),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure"))
{
    for (size_t j = 0; j != contact_particles_.size(); ++j)
    {
        dv_contact_Vol_.push_back(
            contact_particles_[j]->template getVariableByName<Real>("VolumetricMeasure"));
        contact_target_indices_.push_back(j);
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class TargetIdentifier>
Interaction<Contact<Parameters...>>::Interaction(
    ContactRelationType &contact_relation, const StdVec<TargetIdentifier *> &target_identifiers)
    : BaseLocalDynamicsType(contact_relation.getSourceIdentifier()),
      contact_relation_(&contact_relation),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure"))
{
    StdVec<SPHBody *> all_contact_bodies = contact_relation.getContactBodies();
    StdVec<BaseParticles *> all_contact_particles = contact_relation.getContactParticles();
    StdVec<SPHAdaptation *> all_contact_adaptations = contact_relation.getContactAdaptations();
    for (size_t k = 0; k != target_identifiers.size(); ++k)
    {
        TargetIdentifier *target_identifier = target_identifiers[k];
        bool target_found = false;
        for (size_t j = 0; j != all_contact_bodies.size(); ++j)
        {
            if (all_contact_bodies[j]->Name() == target_identifier->Name())
            {
                contact_bodies_.push_back(all_contact_bodies[j]);
                contact_particles_.push_back(all_contact_particles[j]);
                contact_adaptations_.push_back(all_contact_adaptations[j]);
                contact_target_indices_.push_back(j);
                dv_contact_Vol_.push_back(
                    all_contact_particles[j]->template getVariableByName<Real>("VolumetricMeasure"));
                target_found = true;
                break;
            }
        }

        if (!target_found)
        {
            std::cout << "Error: target body " << target_identifier->Name()
                      << " is not found in contact relation " << contact_relation.Name() << std::endl;
            std::exit(1);
        }
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class TargetIdentifier>
Interaction<Contact<Parameters...>>::Interaction(
    ContactRelationType &contact_relation, TargetIdentifier &target_identifiers)
    : Interaction(contact_relation, StdVec<TargetIdentifier *>({&target_identifiers})) {}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Contact<Parameters...>>::
    registerComputingKernel(Implementation<Base> *implementation, UnsignedInt target_index)
{
    contact_relation_->registerComputingKernel(implementation, getContactIndex(target_index));
}
//=================================================================================================//
template <typename... Parameters>
void Interaction<Contact<Parameters...>>::resetComputingKernelUpdated(UnsignedInt target_index)
{
    contact_relation_->resetComputingKernelUpdated(getContactIndex(target_index));
}
//=================================================================================================//
template <typename... Parameters>
UnsignedInt Interaction<Contact<Parameters...>>::getContactIndex(UnsignedInt target_index)
{
    return contact_target_indices_[target_index];
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
Interaction<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt target_index)
    : NeighborList(ex_policy, *encloser.contact_relation_, encloser.getContactIndex(target_index)),
      NeighborKernel(ex_policy, encloser.contact_relation_
                                    ->getNeighborhood(encloser.getContactIndex(target_index))) {}
//=================================================================================================//
template <class WallContactRelationType>
Interaction<Wall>::Interaction(WallContactRelationType &wall_contact_relation)
{
    StdVec<BaseParticles *> contact_particles = wall_contact_relation.getContactParticles();
    StdVec<SPHBody *> contact_bodies = wall_contact_relation.getContactBodies();
    for (size_t k = 0; k != contact_particles.size(); ++k)
    {
        Solid &solid_material = DynamicCast<Solid>(this, contact_bodies[k]->getMatterMaterial());
        dv_wall_vel_ave_.push_back(solid_material.AverageVelocityVariable(contact_particles[k]));
        dv_wall_acc_ave_.push_back(solid_material.AverageAccelerationVariable(contact_particles[k]));
        dv_wall_n_.push_back(contact_particles[k]->template getVariableByName<Vecd>("NormalDirection"));
    }
}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_CK_HPP
