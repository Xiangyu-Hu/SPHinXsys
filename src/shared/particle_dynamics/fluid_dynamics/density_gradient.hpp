#pragma once

#include "density_gradient.h"

namespace SPH
{
namespace fluid_dynamics
{

template <class DataDelegationType>
template <class BaseRelationType>
DensityGradient<DataDelegationType>::DensityGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      rho_(this->particles_->template getVariableDataByName<Real>("Density"))
{
    this->particles_->template registerStateVariable<Vecd>("DensityGradient");
    rho_grad_ = this->particles_->template getVariableDataByName<Vecd>("DensityGradient");
}

template <class KernelCorrectionType>
DensityGradient<Inner<KernelCorrectionType>>::DensityGradient(BaseInnerRelation &inner_relation)
    : DensityGradient<DataDelegateInner>(inner_relation),
      kernel_correction_(particles_) {}

template <class KernelCorrectionType>
void DensityGradient<Inner<KernelCorrectionType>>::interaction(size_t index_i, Real dt)
{
    Vecd density_grad = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real rho_i = rho_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real rho_j = rho_[index_j];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        density_grad -= (rho_i - rho_j) * nablaW_ijV_j;
    }

    rho_grad_[index_i] = density_grad;
}

template <class KernelCorrectionType>
void DensityGradient<Inner<KernelCorrectionType>>::update(size_t index_i, Real dt)
{
    rho_grad_[index_i] = kernel_correction_(index_i) * rho_grad_[index_i];
}

} // namespace fluid_dynamics
} // namespace SPH