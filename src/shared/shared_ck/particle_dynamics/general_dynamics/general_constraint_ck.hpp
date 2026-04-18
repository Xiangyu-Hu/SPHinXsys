#ifndef GENERAL_CONSTRAINT_CK_HPP
#define GENERAL_CONSTRAINT_CK_HPP

#include "general_constraint_ck.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsIdentifier>
MotionConstraintCK<DynamicsIdentifier>::MotionConstraintCK(DynamicsIdentifier &identifier)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_pos0_(this->particles_->template registerStateVariableFrom<Vecd>("InitialPosition", "Position")),
      dv_vel_(this->particles_->template registerStateVariable<Vecd>("Velocity")) {}
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
FixConstraintCK<DynamicsIdentifier>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier>
void FixConstraintCK<DynamicsIdentifier>::UpdateKernel::update(size_t index_i, Real dt)
{
    pos_[index_i] = pos0_[index_i];
    vel_[index_i] = Vecd::Zero();
};
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_CONSTRAINT_CK_HPP
