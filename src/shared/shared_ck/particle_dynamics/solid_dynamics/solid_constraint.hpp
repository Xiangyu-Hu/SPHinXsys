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
      dv_n_(this->particles_->template getVariableByName<Vecd>("NormalDirection")),
      dv_n0_(this->particles_->template registerStateVariableOnlyFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      dv_acc_(this->particles_->template registerStateVariableOnly<Vecd>("Acceleration")),
      sv_simbody_state_(this->particles_->template addUniqueSingularVariableOnly<SimbodyState>("SimbodyState"))
{
    this->particles_->template addVariableToWrite<Vecd>("Velocity");
    const SimTK::State *state = &integ_.getState();
    MBsystem_.realize(*state, SimTK::Stage::Acceleration);
    sim_tk_initial_origin_location_ = mobod_.getBodyOriginLocation(*state);
    sv_simbody_state_->setValue(SimbodyState(sim_tk_initial_origin_location_, mobod_, *state));
}
//=================================================================================================//
template <class DynamicsIdentifier>
void ConstraintBySimBodyCK<DynamicsIdentifier>::setupDynamics(Real dt)
{
    const SimTK::State *state = &integ_.getState();
    MBsystem_.realize(*state, SimTK::Stage::Acceleration);
    sv_simbody_state_->setValue(SimbodyState(sim_tk_initial_origin_location_, mobod_, *state));
};
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
ConstraintBySimBodyCK<DynamicsIdentifier>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      pos0_(encloser.dv_pos0_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      n0_(encloser.dv_n0_->DelegatedData(ex_policy)),
      acc_(encloser.dv_acc_->DelegatedData(ex_policy)),
      simbody_state_(encloser.sv_simbody_state_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier>
void ConstraintBySimBodyCK<DynamicsIdentifier>::UpdateKernel::update(size_t index_i, Real dt)
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
template <class DynamicsIdentifier>
TotalForceForSimBodyCK<DynamicsIdentifier>::
    TotalForceForSimBodyCK(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                           SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ)
    : BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>(identifier),
      MBsystem_(MBsystem), mobod_(mobod), integ_(integ),
      dv_force_(this->particles_->template registerStateVariableOnly<Vecd>("Force")),
      dv_force_prior_(this->particles_->template getVariableByName<Vecd>("ForcePrior")),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      sv_current_origin_location_(
          this->particles_->template addUniqueSingularVariableOnly<Vec3d>(
              identifier.getName() + "OriginLocation"))
{
    this->quantity_name_ = "TotalForceForSimBody";
}
//=================================================================================================//
template <class DynamicsIdentifier>
void TotalForceForSimBodyCK<DynamicsIdentifier>::setupDynamics(Real dt)
{
    const SimTK::State *state = &integ_.getState();
    MBsystem_.realize(*state, SimTK::Stage::Acceleration);
    sv_current_origin_location_->setValue(SimTKToEigen(mobod_.getBodyOriginLocation(*state)));
}
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
TotalForceForSimBodyCK<DynamicsIdentifier>::
    ReduceKernel::ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      current_origin_location_(encloser.sv_current_origin_location_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier>
SimTK::SpatialVec TotalForceForSimBodyCK<DynamicsIdentifier>::
    ReduceKernel::reduce(size_t index_i, Real dt)
{
    Vecd force = force_[index_i] + force_prior_[index_i];
    Vec3d force_from_particle = upgradeToVec3d(force);
    Vecd displacement = pos_[index_i] - degradeToVecd(*current_origin_location_);
    Vec3d torque_from_particle = upgradeToVec3d(displacement).cross(force_from_particle);
    return SimTK::SpatialVec(EigenToSimTK(torque_from_particle), EigenToSimTK(force_from_particle));
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // SOLID_CONSTRAINT_HPP
