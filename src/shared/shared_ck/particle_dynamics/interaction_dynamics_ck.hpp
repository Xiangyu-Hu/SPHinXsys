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
      pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      kernel_(encloser.kernel_) {}
//=================================================================================================//
template <class ExecutionPolicy>
InteractionDynamics<Contact<>>::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    InteractionDynamics<Inner<>> &encloser, UnsignedInt contact_index)
    : NeighborList(ex_policy,
                   encloser.dv_contact_neighbor_index_[contact_index],
                   encloser.dv_contact_particle_offset_[contact_index]),
      pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      contact_pos_(encloser.dv_contact_pos_[contact_index]->DelegatedDataField(ex_policy)),
      kernel_(encloser.contact_kernel_[contact_index]) {}
//=================================================================================================//
} // namespace SPH
#endif // INTERACTION_DYNAMICS_CK_HPP