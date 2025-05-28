#pragma once

#include "energy_gradient.h"

namespace SPH
{
namespace fluid_dynamics
{

template <class DataDelegationType>
template <class BaseRelationType>
EnergyGradient<DataDelegationType>::EnergyGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      energy_(this->particles_->template getVariableDataByName<Real>("Energy")),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      energy_grad_(this->particles_->template registerStateVariable<Vecd>("EnergyGradient")) {}
template <class KernelCorrectionType>
EnergyGradient<Inner<KernelCorrectionType>>::EnergyGradient(BaseInnerRelation &inner_relation)
    : EnergyGradient<DataDelegateInner>(inner_relation),
      kernel_correction_(particles_) {}

template <class KernelCorrectionType>
void EnergyGradient<Inner<KernelCorrectionType>>::interaction(size_t index_i, Real dt)
{
    Vecd energy_grad = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real energy_i = energy_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real energy_j = energy_[index_j];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        energy_grad -= (energy_i - energy_j) * nablaW_ijV_j;
    }

    energy_grad_[index_i] = energy_grad;
}

template <class KernelCorrectionType>
void EnergyGradient<Inner<KernelCorrectionType>>::update(size_t index_i, Real dt)
{
    energy_grad_[index_i] = kernel_correction_(index_i)* energy_grad_[index_i];
}

} // namespace fluid_dynamics
} // namespace SPH
