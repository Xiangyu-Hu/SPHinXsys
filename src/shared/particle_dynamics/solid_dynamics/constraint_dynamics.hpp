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
      physical_time_(this->sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")),
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
      n0_(this->particles_->template registerSharedVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      acc_(this->particles_->template registerSharedVariable<Vecd>("Acceleration"))
{
    simbody_state_ = &integ_.getState();
    MBsystem_.realize(*simbody_state_, SimTK::Stage::Acceleration);
    initial_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
}
//=================================================================================================//
template <class DynamicsIdentifier>
void ConstraintBySimBody<DynamicsIdentifier>::update(size_t index_i, Real dt)
{
    /** Change to SimTK::Vector. */
    SimTKVec3 rr, pos, vel, acc;
    rr = EigenToSimTK(upgradeToVec3d(this->pos0_[index_i])) - initial_mobod_origin_location_;
    mobod_.findStationLocationVelocityAndAccelerationInGround(*simbody_state_, rr, pos, vel, acc);
    /** this is how we calculate the particle position in after transform of MBbody.
     * const SimTK::Rotation&  R_GB = mobod_.getBodyRotation(simbody_state);
     * const SimTKVec3&      p_GB = mobod_.getBodyOriginLocation(simbody_state);
     * const SimTKVec3 r = R_GB * rr; // re-express station vector p_BS in G (15 flops)
     * base_particle_data_i.pos_  = (p_GB + r);
     */
    this->pos_[index_i] = degradeToVecd(SimTKToEigen(pos));
    this->vel_[index_i] = degradeToVecd(SimTKToEigen(vel));
    acc_[index_i] = degradeToVecd(SimTKToEigen(acc));

    SimTKVec3 n = (mobod_.getBodyRotation(*simbody_state_) * EigenToSimTK(upgradeToVec3d(n0_[index_i])));
    n_[index_i] = degradeToVecd(SimTKToEigen(n));
}
//=================================================================================================//
template <class DynamicsIdentifier>
TotalForceForSimBody<DynamicsIdentifier>::
    TotalForceForSimBody(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                         SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ)
    : BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>(identifier),
      DataDelegateSimple(identifier.getSPHBody()),
      force_(particles_->registerSharedVariable<Vecd>("Force")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      MBsystem_(MBsystem), mobod_(mobod), integ_(integ)
{
    this->quantity_name_ = "TotalForceForSimBody";
}
//=================================================================================================//
template <class DynamicsIdentifier>
void TotalForceForSimBody<DynamicsIdentifier>::setupDynamics(Real dt)
{
    const SimTK::State *simbody_state = &integ_.getState();
    MBsystem_.realize(*simbody_state, SimTK::Stage::Acceleration);
    current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state);
}
//=================================================================================================//
template <class DynamicsIdentifier>
SimTK::SpatialVec TotalForceForSimBody<DynamicsIdentifier>::reduce(size_t index_i, Real dt = 0.0)
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
