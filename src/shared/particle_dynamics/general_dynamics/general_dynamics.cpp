#include "general_dynamics.h"
#include <limits>

namespace SPH
{
//=================================================================================================//
TimeStepInitialization::TimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : BaseTimeStepInitialization(sph_body, gravity_ptr), GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_), acc_prior_(particles_->acc_prior_) {}
//=================================================================================================//
void TimeStepInitialization::update(size_t index_i, Real dt)
{
    acc_prior_[index_i] = gravity_->InducedAcceleration(pos_[index_i]);
}
//=================================================================================================//
RandomizeParticlePosition::RandomizeParticlePosition(SPHBody &sph_body)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_), randomize_scale_(sph_body.sph_adaptation_->MinimumSpacing()) {}
//=================================================================================================//
void RandomizeParticlePosition::update(size_t index_i, Real dt)
{
    Vecd &pos_n_i = pos_[index_i];
    for (int k = 0; k < pos_n_i.size(); ++k)
    {
        pos_n_i[k] += dt * (((Real)rand() / (RAND_MAX)) - 0.5) * 2.0 * randomize_scale_;
    }
}
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
UpperFrontInXDirection::UpperFrontInXDirection(SPHBody &sph_body)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_)
{
    quantity_name_ = "UpperFrontInXDirection";
}
//=================================================================================================//
Real UpperFrontInXDirection::reduce(size_t index_i, Real dt)
{
    return pos_[index_i][0];
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
    : LocalDynamicsReduce<Vecd, ReduceLowerBound>(sph_body, MaxRealNumber * Vecd::Ones()),
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
    : LocalDynamicsReduce<Vecd, ReduceUpperBound>(sph_body, MinRealNumber * Vecd::Ones()),
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
TotalMechanicalEnergy::TotalMechanicalEnergy(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : LocalDynamicsReduce<Real, ReduceSum<Real>>(sph_body, Real(0)),
      GeneralDataDelegateSimple(sph_body), mass_(particles_->mass_),
      vel_(particles_->vel_), pos_(particles_->pos_),
      gravity_(gravity_ptr_keeper_.assignPtr(gravity_ptr))
{
    quantity_name_ = "TotalMechanicalEnergy";
}
//=================================================================================================//
Real TotalMechanicalEnergy::reduce(size_t index_i, Real dt)
{
    return 0.5 * mass_[index_i] * vel_[index_i].squaredNorm() + mass_[index_i] * gravity_->getPotential(pos_[index_i]);
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//