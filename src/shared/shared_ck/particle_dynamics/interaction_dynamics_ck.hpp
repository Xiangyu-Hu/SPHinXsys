#ifndef INTERACTION_DYNAMICS_CK_HPP
#define INTERACTION_DYNAMICS_CK_HPP

#include "interaction_dynamics_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
InteractionDynamics<Inner<>>::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy, InteractionDynamics<Inner<>> &encloser)
    : NeighborList(ex_policy, encloser.dv_neighbor_index_, encloser.dv_particle_offset_),
      neighbor_(ex_policy, encloser.sph_adaptation_, encloser.dv_pos_) {}
//=================================================================================================//
template <class ExecutionPolicy>
InteractionDynamics<Contact<>>::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    InteractionDynamics<Contact<>> &encloser, UnsignedInt contact_index)
    : NeighborList(ex_policy,
                   encloser.dv_contact_neighbor_index_[contact_index],
                   encloser.dv_contact_particle_offset_[contact_index]),
      neighbor_(ex_policy, encloser.sph_adaptation_, encloser.contact_adaptations_[contact_index],
                encloser.dv_pos_, encloser.contact_pos_[contact_index]) {}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_DYNAMICS_CK_HPP