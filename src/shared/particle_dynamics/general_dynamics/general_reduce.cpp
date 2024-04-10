#include "general_reduce.h"
#include <limits>

namespace SPH
{
//=================================================================================================//
VelocityBoundCheck::
    VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound)
    : LocalDynamicsReduce<bool, ReduceOR>(sph_body, false),
      GeneralDataDelegateSimple(sph_body),
      vel_(particles_->vel_), velocity_bound_(velocity_bound) {}
//=================================================================================================//
bool VelocityBoundCheck::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].norm() > velocity_bound_;
}
//=================================================================================================//
MaximumSpeed::MaximumSpeed(SPHBody &sph_body)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      GeneralDataDelegateSimple(sph_body),
      vel_(particles_->vel_)
{
    quantity_name_ = "MaximumSpeed";
}
//=================================================================================================//
Real MaximumSpeed::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].norm();
}
//=================================================================================================//
PositionLowerBound::PositionLowerBound(SPHBody &sph_body)
    : LocalDynamicsReduce<Vecd, ReduceLowerBound>(sph_body, MaxReal * Vecd::Ones()),
      GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_)
{
    quantity_name_ = "PositionLowerBound";
}
//=================================================================================================//
Vecd PositionLowerBound::reduce(size_t index_i, Real dt)
{
    return pos_[index_i];
}
//=================================================================================================//
PositionUpperBound::PositionUpperBound(SPHBody &sph_body)
    : LocalDynamicsReduce<Vecd, ReduceUpperBound>(sph_body, MinReal * Vecd::Ones()),
      GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_)
{
    quantity_name_ = "PositionUpperBound";
}
//=================================================================================================//
Vecd PositionUpperBound::reduce(size_t index_i, Real dt)
{
    return pos_[index_i];
}
//=================================================================================================//
TotalKineticEnergy::TotalKineticEnergy(SPHBody &sph_body)
    : LocalDynamicsReduce<Real, ReduceSum<Real>>(sph_body, Real(0)),
      GeneralDataDelegateSimple(sph_body),
      mass_(particles_->mass_), vel_(particles_->vel_)
{
    quantity_name_ = "TotalKineticEnergy";
}
//=================================================================================================//
Real TotalKineticEnergy::reduce(size_t index_i, Real dt)
{
    return 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
}
//=================================================================================================//
TotalMechanicalEnergy::TotalMechanicalEnergy(SPHBody &sph_body, Gravity &gravity)
    : TotalKineticEnergy(sph_body),
      gravity_(gravity), pos_(particles_->pos_),
      device_kernel(particles_, gravity_)
{
    quantity_name_ = "TotalMechanicalEnergy";
}
//=================================================================================================//
Real TotalMechanicalEnergy::reduce(size_t index_i, Real dt)
{
    return TotalKineticEnergy::reduce(index_i, dt) +
           mass_[index_i] * gravity_.getPotential(pos_[index_i]);
}
//=================================================================================================//
} // namespace SPH
