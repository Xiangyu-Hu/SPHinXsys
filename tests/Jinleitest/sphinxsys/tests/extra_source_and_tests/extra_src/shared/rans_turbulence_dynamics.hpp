#ifndef RANS_TURBULENCE_DYNAMICS_HPP
#define RANS_TURBULENCE_DYNAMICS_HPP
#include "rans_turbulence_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TurbulentViscousForceInFVM<DataDelegationType>::TurbulentViscousForceInFVM(BaseRelationType &base_relation)
    : ForcePrior(base_relation.getSPHBody(), "TurbulentViscousForce"), DataDelegationType(base_relation),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      mu_t_(this->particles_->template getVariableDataByName<Real>("TurblunetViscosity")),
      wall_adjacent_cell_flag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      turbulent_viscous_force_(this->particles_->template registerStateVariable<Vecd>("TurbulentViscousForce")),
      smoothing_length_(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
TurbulentViscousForceInFVM<Inner<>>::TurbulentViscousForceInFVM(BaseInnerRelation &inner_relation)
    : TurbulentViscousForceInFVM<DataDelegateInner>(inner_relation)
{
}
//=================================================================================================//
void TurbulentViscousForceInFVM<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

        // turbulent viscous force
        Vecd vel_derivative = (vel_[index_i] - vel_[index_j]) /
                              (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);

        force += e_ij.dot(2.0 * e_ij) * mu_t_avg * vel_derivative * inner_neighborhood.dW_ij_[n] * Vol_[index_j];
    }
    turbulent_viscous_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
TkeGradientForceInFVM<DataDelegationType>::TkeGradientForceInFVM(BaseRelationType &base_relation)
    : ForcePrior(base_relation.getSPHBody(), "TkeGradientForce"), DataDelegationType(base_relation),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      K_(this->particles_->template getVariableDataByName<Real>("TKE")),
      tke_gradient_force_(this->particles_->template registerStateVariable<Vecd>("TkeGradientForce"))
{
}
//=================================================================================================//
TkeGradientForceInFVM<Inner<>>::TkeGradientForceInFVM(BaseInnerRelation &inner_relation)
    : TkeGradientForceInFVM<DataDelegateInner>(inner_relation)
{
}
//=================================================================================================//
void TkeGradientForceInFVM<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];

        // tke gradient force
        force += (inner_neighborhood.dW_ij_[n] * Vol_[index_j] * rho_[index_i] * (2.0 / 3.0) * (K_[index_i] - K_[index_j]) * Matd::Identity()) * e_ij;
    }

    tke_gradient_force_[index_i] = force * Vol_[index_i];
}
} // namespace fluid_dynamics
} // namespace SPH
#endif // RANS_TURBULENCE_DYNAMICS_HPP
