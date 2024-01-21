
#include "eulerian_compressible_fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
EulerianCompressibleTimeStepInitialization::
    EulerianCompressibleTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : TimeStepInitialization(sph_body, gravity_ptr), rho_(particles_->rho_), mass_(particles_->mass_),
      pos_(particles_->pos_), vel_(particles_->vel_),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")),
      dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")){};
//=================================================================================================//
void EulerianCompressibleTimeStepInitialization::update(size_t index_i, Real dt)
{
    dmom_dt_prior_[index_i] = mass_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
    dE_dt_prior_[index_i] = mass_[index_i] * (gravity_->InducedAcceleration(pos_[index_i])).dot(vel_[index_i]);
}
//=================================================================================================//
EulerianCompressibleAcousticTimeStepSize::
    EulerianCompressibleAcousticTimeStepSize(SPHBody &sph_body)
    : AcousticTimeStepSize(sph_body), rho_(particles_->rho_),
      p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_),
      smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()),
      compressible_fluid_(CompressibleFluid(1.0, 1.4)){};
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return 0.6 / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
}
//=================================================================================================//
EulerianCompressibleViscousForceInner::
    EulerianCompressibleViscousForceInner(BaseInnerRelation &inner_relation)
    : ViscousForceInner(inner_relation),
      dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
void EulerianCompressibleViscousForceInner::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    const Vecd &vel_i = vel_[index_i];

    Vecd force = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        // viscous force
        vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        force += 2.0 * mass_[index_i] * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
    }
    dmom_dt_prior_[index_i] += force;
    dE_dt_prior_[index_i] += force.dot(vel_[index_i]);
}
//=================================================================================================//
BaseIntegrationInCompressible::BaseIntegrationInCompressible(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation),
      compressible_fluid_(CompressibleFluid(1.0, 1.4)),
      Vol_(particles_->Vol_), E_(*particles_->getVariableByName<Real>("TotalEnergy")),
      dE_dt_(*particles_->getVariableByName<Real>("TotalEnergyChangeRate")),
      dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")),
      dmass_dt_(*this->particles_->template registerSharedVariable<Real>("MassChangeRate")),
      mom_(*particles_->getVariableByName<Vecd>("Momentum")),
      dmom_dt_(*particles_->getVariableByName<Vecd>("MomentumChangeRate")),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
