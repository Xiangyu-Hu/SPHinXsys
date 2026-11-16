#pragma once

#include "udf_common_turbulence_model.h"

namespace SPH
{

namespace fluid_dynamics
{

namespace udf
{

template <class DataDelegationType>
template <class BaseRelationType>
TKEnergyForce<Base, DataDelegationType>::
    TKEnergyForce(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      force_(this->particles_->template registerStateVariableData<Vecd>("Force")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      indicator_(this->particles_->template getVariableDataByName<int>("Indicator")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      turbu_k_(this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      test_k_grad_rslt_(this->particles_->template registerStateVariableData<Vecd>("TkeGradResult")) {}

template <class DataDelegationType>
template <class BaseRelationType>
kEpsilon_GetVelocityGradient<DataDelegationType>::
    kEpsilon_GetVelocityGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      is_near_wall_P1_(this->particles_->template getVariableDataByName<int>("IsNearWallP1")),
      is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      velocity_gradient_(this->particles_->template registerStateVariableData<Matd>("TurbulentVelocityGradient")),
      velocity_gradient_wall(this->particles_->template registerStateVariableData<Matd>("Velocity_Gradient_Wall")) {}

template <class DataDelegationType>
template <class BaseRelationType>
TurbuViscousForce<DataDelegationType>::TurbuViscousForce(BaseRelationType &base_relation)
    : ViscousForce<DataDelegationType>(base_relation),

      turbu_k_(this->particles_->template getVariableDataByName<Real>("TurbulenceKineticEnergy")),
      turbu_mu_(this->particles_->template getVariableDataByName<Real>("TurbulentViscosity")),
      wall_Y_plus_(this->particles_->template getVariableDataByName<Real>("WallYplus")),
      wall_Y_star_(this->particles_->template getVariableDataByName<Real>("WallYstar")),
      velo_friction_(this->particles_->template getVariableDataByName<Vecd>("FrictionVelocity")),
      y_p_(this->particles_->template getVariableDataByName<Real>("Y_P")),
      is_near_wall_P2_(this->particles_->template getVariableDataByName<int>("IsNearWallP2")),
      viscosity_(this->sph_body_-> template getMaterialProperty<Viscosity>()),
      molecular_viscosity_(viscosity_.ReferenceViscosity()),
      c0_(DynamicCast<Fluid>(this, this->sph_body_->getMatterMaterial()).ReferenceSoundSpeed()) {}

template <class DataDelegationType>
template <class BaseRelationType>
TurbulentLinearGradientCorrectionMatrix<DataDelegationType>::
    TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      turbu_B_(this->particles_->template registerStateVariableData<Matd>(
          "TurbulentLinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)),
      B_only_wall_(this->particles_->template registerStateVariableData<Matd>(
          "TurbulentLinearGradientCorrectionMatrixOnlyWall", IdentityMatrix<Matd>::value))
{

}

template <class RiemannSolverType>
TurbulentIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::
    TurbulentIntegration2ndHalf(BaseContactRelation &wall_contact_relation)
    : BaseIntegrationWithWall(wall_contact_relation),
      riemann_solver_(this->fluid_, this->fluid_, 3.0) {}

template <class RiemannSolverType>
void TurbulentIntegration2ndHalf<Contact<Wall>, RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd *vel_ave_k = wall_vel_ave_[k];
        Vecd *n_k = wall_n_[k];
        Real *wall_Vol_k = wall_Vol_[k];
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];

            Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - vel_[index_i];
            density_change_rate += (vel_[index_i] - vel_in_wall).dot(e_ij) * dW_ijV_j;
            Real u_jump = 2.0 * (vel_[index_i] - vel_ave_k[index_j]).dot(n_k[index_j]);
            p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * n_k[index_j];
        }
    }
    drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    force_[index_i] += p_dissipation * this->Vol_[index_i];
}

}

}

}
