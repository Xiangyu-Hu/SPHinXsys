#pragma once
#include "shear_dynamics.h"
namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseShearStressIntegration<DataDelegationType>::BaseShearStressIntegration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      vel_(this->particles_->vel_),
      velocity_gradient_(this->particles_->velocity_gradient_),
      shear_stress_(this->particles_->shear_stress_),
      p_(*this->particles_->template getVariableByName<Real>("Pressure")),
      von_mises_stress_(this->particles_->von_mises_stress_) {}
//=================================================================================================//
template <class DataDelegationType>
void BaseShearStressIntegration<DataDelegationType>::interaction(size_t index_i, Real dt)
{
    Matd velocity_gradient = Matd::Zero();
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        Matd velocity_gradient_ij = -(vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
        velocity_gradient += velocity_gradient_ij;
    }
    velocity_gradient_[index_i] = velocity_gradient;
}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseShearStressAcceleration<DataDelegationType>::BaseShearStressAcceleration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      shear_stress_(*this->particles_->template getVariableByName<Matd>("ShearStress")),
      rho_(this->particles_->rho_), acc_prior_(this->particles_->acc_prior_) {}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH