#ifndef CONSTRAINT_DYNAMICS_HPP
#define CONSTRAINT_DYNAMICS_HPP

#include "constraint_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <class DynamicsIdentifier>
PositionTranslate<DynamicsIdentifier>::
    PositionTranslate(DynamicsIdentifier &identifier, Real start_time, Real end_time, Vecd translation)
    : MotionConstraint<DynamicsIdentifier>(identifier),
      start_time_(start_time), end_time_(end_time),
      physical_time_(this->sph_system_.template getSystemVariableDataByName<Real>("PhysicalTime")),
      translation_(translation) {}
//=================================================================================================//
template <class DynamicsIdentifier>
void PositionTranslate<DynamicsIdentifier>::update(size_t index_i, Real dt)
{
    // only apply in the defined time period
    if (*physical_time_ >= start_time_ && *physical_time_ <= end_time_)
    {
        // displacement from the initial position, 0.5x because it's executed twice
        this->pos_[index_i] += 0.5 * getDisplacement(index_i, dt);
        this->vel_[index_i] = Vecd::Zero();
    }
}
//=================================================================================================//
template <class DynamicsIdentifier>
Vecd PositionTranslate<DynamicsIdentifier>::getDisplacement(size_t index_i, Real dt)
{
    return (this->pos0_[index_i] + translation_ - this->pos_[index_i]) * dt /
           (end_time_ - *physical_time_);
}
//=================================================================================================//
template <class DynamicsIdentifier>
ConstraintBySimBody<DynamicsIdentifier>::
    ConstraintBySimBody(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                        SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ)
    : MotionConstraint<DynamicsIdentifier>(identifier),
      MBsystem_(MBsystem), mobod_(mobod), integ_(integ),
      n_(this->particles_->template getVariableDataByName<Vecd>("NormalDirection")),
      n0_(this->particles_->template registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      acc_(this->particles_->template registerStateVariable<Vecd>("Acceleration"))
{
    const SimTK::State *state = &integ_.getState();
    MBsystem_.realize(*state, SimTK::Stage::Acceleration);
    initializeSimbodyState(*state);
}
//=================================================================================================//
template <class DynamicsIdentifier>
void ConstraintBySimBody<DynamicsIdentifier>::initializeSimbodyState(const SimTK::State &state)
{
    updateSimbodyState(state);
    simbody_state_.initial_origin_location_ = simbody_state_.origin_location_;
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
void ConstraintBySimBody<DynamicsIdentifier>::updateSimbodyState(const SimTK::State &state)
{
    simbody_state_.origin_location_ = SimTKToEigen(mobod_.getBodyOriginLocation(state));
    simbody_state_.origin_velocity_ = SimTKToEigen(mobod_.getBodyOriginVelocity(state));
    simbody_state_.origin_acceleration_ = SimTKToEigen(mobod_.getBodyOriginAcceleration(state));
    simbody_state_.angular_velocity_ = SimTKToEigen(mobod_.getBodyAngularVelocity(state));
    simbody_state_.angular_acceleration_ = SimTKToEigen(mobod_.getBodyAngularAcceleration(state));
    simbody_state_.rotation_ = SimTKToEigen(mobod_.getBodyRotation(state));
}
//=================================================================================================//
template <class DynamicsIdentifier>
void ConstraintBySimBody<DynamicsIdentifier>::update(size_t index_i, Real dt)
{
    Vec3d pos, vel, acc, n;
    simbody_state_.findStationLocationVelocityAndAccelerationInGround(
        upgradeToVec3d(this->pos0_[index_i]), upgradeToVec3d(n0_[index_i]), pos, vel, acc, n);
    this->pos_[index_i] = degradeToVecd(pos);
    this->vel_[index_i] = degradeToVecd(vel);
    acc_[index_i] = degradeToVecd(acc);
    n_[index_i] = degradeToVecd(n);
}
//=================================================================================================//
template <class DynamicsIdentifier>
TotalForceForSimBody<DynamicsIdentifier>::
    TotalForceForSimBody(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                         SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ)
    : BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>(identifier),

      force_(this->particles_->template registerStateVariable<Vecd>("Force")),
      force_prior_(this->particles_->template getVariableDataByName<Vecd>("ForcePrior")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      MBsystem_(MBsystem), mobod_(mobod), integ_(integ)
{
    this->quantity_name_ = "TotalForceForSimBody";
}
//=================================================================================================//
template <class DynamicsIdentifier>
void TotalForceForSimBody<DynamicsIdentifier>::setupDynamics(Real dt)
{
    const SimTK::State *state = &integ_.getState();
    MBsystem_.realize(*state, SimTK::Stage::Acceleration);
    current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*state);
}
//=================================================================================================//
template <class DynamicsIdentifier>
SimTK::SpatialVec TotalForceForSimBody<DynamicsIdentifier>::reduce(size_t index_i, Real dt)
{
    Vecd force = force_[index_i] + force_prior_[index_i];
    SimTKVec3 force_from_particle = EigenToSimTK(upgradeToVec3d(force));
    SimTKVec3 displacement = EigenToSimTK(upgradeToVec3d(pos_[index_i])) - current_mobod_origin_location_;
    SimTKVec3 torque_from_particle = SimTK::cross(displacement, force_from_particle);
    return SimTK::SpatialVec(torque_from_particle, force_from_particle);
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // CONSTRAINT_DYNAMICS_HPP
