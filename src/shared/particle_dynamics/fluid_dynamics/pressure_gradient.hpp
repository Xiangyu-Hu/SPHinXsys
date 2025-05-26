// pressure_gradient.hpp
#pragma once

#include "pressure_gradient.h"

namespace SPH
{
namespace fluid_dynamics
{

template <class DataDelegationType>
template <class BaseRelationType>
PressureGradient<DataDelegationType>::PressureGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      p_(this->particles_->template getVariableDataByName<Real>("Pressure")),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")), // YOĞUNLUK EKLENDİ
      p_grad_(this->particles_->template registerStateVariable<Vecd>("PressureGradient")) {}

template <class KernelCorrectionType>
PressureGradient<Inner<KernelCorrectionType>>::PressureGradient(BaseInnerRelation &inner_relation)
    : PressureGradient<DataDelegateInner>(inner_relation),
      kernel_correction_(particles_) {}

template <class KernelCorrectionType>
void PressureGradient<Inner<KernelCorrectionType>>::interaction(size_t index_i, Real dt)
{
    Vecd pressure_grad = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real p_i = p_[index_i];

    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real p_j = p_[index_j];

        // Follow exactly the same pattern as velocity gradient
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        
        // Note: For pressure (scalar), we don't need transpose() as in velocity (vector)
        pressure_grad -= (p_i - p_j) * nablaW_ijV_j;
    }

    // Direct assignment as in velocity gradient
    p_grad_[index_i] = pressure_grad;
}

template <class KernelCorrectionType>
void PressureGradient<Inner<KernelCorrectionType>>::update(size_t index_i, Real dt)
{
    // Get the kernel correction matrix
    Matd correction_matrix = kernel_correction_(index_i);
    // Extract the diagonal elements and use them as correction factors
    Vecd correction_factors;
    for (int i = 0; i < Dimensions; ++i)
        correction_factors[i] = correction_matrix(i,i);
    
    // Apply correction component-wise
    p_grad_[index_i] = correction_factors.cwiseProduct(p_grad_[index_i]);
}

} // namespace fluid_dynamics
} // namespace SPH
