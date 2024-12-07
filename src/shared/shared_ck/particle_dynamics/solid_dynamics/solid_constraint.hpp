#ifndef SOLID_CONSTRAINT_HPP
#define SOLID_CONSTRAINT_HPP

#include "solid_constraint.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <class DynamicsIdentifier>
ConstraintBySimBodyCK<DynamicsIdentifier>::
    ConstraintBySimBodyCK(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                          SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      MBsystem_(MBsystem), mobod_(mobod), integ_(integ),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_pos0_(this->particles_->template registerStateVariableOnlyFrom<Vecd>("InitialPosition", "Position")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      dv_n_(this->particles_->template getVariableByName<Vecd>("Normal")),
      dv_n0_(this->particles_->template registerStateVariableOnlyFrom<Vecd>("InitialNormal", "Normal")),
      dv_acc_(this->particles_->template getVariableByName<Vecd>("Acceleration")),
      sv_simbody_state_(this->particles_->template addUniqueSingularVariableOnly<SimbodyState>("SimbodyState"))
{
    const SimTK::State *state = &integ_.getState();
    MBsystem_.realize(*state, SimTK::Stage::Acceleration);
    initializeSimbodyState(*state);
}
//=================================================================================================//
template <class DynamicsIdentifier>
void ConstraintBySimBody<DynamicsIdentifier>::setupDynamics(Real dt)
{
    const SimTK::State *state = &integ_.getState();
    MBsystem_.realize(*state, SimTK::Stage::Acceleration);
    updateSimbodyState(*state);
};
//=================================================================================================//
template <class DynamicsIdentifier>
void ConstraintBySimBodyCK<DynamicsIdentifier>::initializeSimbodyState(const SimTK::State &state)
{
    updateSimbodyState(state);
    SimbodyState *simbody_state = sv_simbody_state_->ValueAddress();
    simbody_state->initial_origin_location_ = simbody_state->origin_location_;
}
//=================================================================================================//
void ConstraintBySimBodyCK<DynamicsIdentifier>::updateSimbodyState(const SimTK::State &state)
{
    SimbodyState *simbody_state = sv_simbody_state_->ValueAddress();
    simbody_state->origin_location_ = SimTKToEigen(mobod_.getBodyOriginLocation(state));
    simbody_state->origin_velocity_ = SimTKToEigen(mobod_.getBodyOriginVelocity(state));
    simbody_state->origin_acceleration_ = SimTKToEigen(mobod_.getBodyOriginAcceleration(state));
    simbody_state->angular_velocity_ = SimTKToEigen(mobod_.getBodyAngularVelocity(state));
    simbody_state->angular_acceleration_ = SimTKToEigen(mobod_.getBodyAngularAcceleration(state));
    simbody_state->rotation_ = SimTKToEigen(mobod_.getBodyRotation(state));
}
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
ConstraintBySimBodyCK<DynamicsIdentifier>::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      pos0_(encloser.dv_pos0_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      n_(encloser.dv_n_->DelegatedDataField(ex_policy)),
      n0_(encloser.dv_n0_->DelegatedDataField(ex_policy)),
      acc_(encloser.dv_acc_->DelegatedDataField(ex_policy)),
      simbody_state_(encloser.sv_simbody_state_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier>
ConstraintBySimBodyCK<DynamicsIdentifier>::UpdateKernel::update(size_t index_i, Real dt)
{
    Vec3d pos, vel, acc, n;
    simbody_state_->findStationLocationVelocityAndAccelerationInGround(
        upgradeToVec3d(pos0_[index_i]), upgradeToVec3d(n0_[index_i]), pos, vel, acc, n);
    pos_[index_i] = degradeToVecd(pos);
    vel_[index_i] = degradeToVecd(vel);
    acc_[index_i] = degradeToVecd(acc);
    n_[index_i] = degradeToVecd(n);
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // SOLID_CONSTRAINT_HPP
