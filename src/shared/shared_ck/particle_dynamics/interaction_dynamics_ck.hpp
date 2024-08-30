#ifndef INTERACTION_DYNAMICS_CK_HPP
#define INTERACTION_DYNAMICS_CK_HPP

#include "interaction_dynamics_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename... Parameters>
InteractionDynamics<Inner<Parameters...>>::InteractionDynamics(InnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      sph_adaptation_(sph_body_.sph_adaptation_),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_neighbor_index_(inner_relation.getNeighborIndex()),
      dv_particle_offset_(inner_relation.getParticleOffset()) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
InteractionDynamics<Inner<Parameters...>>::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    InteractionDynamics<Inner<Parameters...>> &encloser)
    : NeighborList(ex_policy, encloser.dv_neighbor_index_, encloser.dv_particle_offset_),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_, encloser.dv_pos_) {}
//=================================================================================================//
template <typename... Parameters>
InteractionDynamics<Contact<Parameters...>>::InteractionDynamics(ContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()),
      sph_adaptation_(sph_body_.sph_adaptation_),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      contact_bodies_(contact_relation.getContactBodies()),
      contact_adaptations_(contact_relation.getContactAdaptations()),
      dv_contact_neighbor_index_(contact_relation.getContactNeighborIndex()),
      dv_contact_particle_offset_(contact_relation.getContactParticleOffset()) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
InteractionDynamics<Contact<Parameters...>>::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    InteractionDynamics<Contact<Parameters...>> &encloser, UnsignedInt contact_index)
    : NeighborList(ex_policy,
                   encloser.dv_contact_neighbor_index_[contact_index],
                   encloser.dv_contact_particle_offset_[contact_index]),
      Neighbor<Parameters...>(ex_policy, encloser.sph_adaptation_,
                              encloser.contact_adaptations_[contact_index],
                              encloser.dv_pos_, encloser.contact_pos_[contact_index]) {}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_DYNAMICS_CK_HPP