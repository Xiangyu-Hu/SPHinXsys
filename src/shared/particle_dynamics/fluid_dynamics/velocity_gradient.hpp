#pragma once

#include "velocity_gradient.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
VelocityGradient<DataDelegationType>::VelocityGradient(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      vel_grad_(this->particles_->template registerStateVariable<Matd>("VelocityGradient")) {}
//=================================================================================================//
template <class KernelCorrectionType>
VelocityGradient<Inner<KernelCorrectionType>>::VelocityGradient(BaseInnerRelation &inner_relation)
    : VelocityGradient<DataDelegateInner>(inner_relation),
      kernel_correction_(particles_) {}
//=================================================================================================//
template <class KernelCorrectionType>
void VelocityGradient<Inner<KernelCorrectionType>>::interaction(size_t index_i, Real dt)
{
    Matd vel_grad = Matd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        vel_grad -= (vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
    }

    vel_grad_[index_i] = vel_grad;
}
//=================================================================================================//
template <class KernelCorrectionType>
void VelocityGradient<Inner<KernelCorrectionType>>::update(size_t index_i, Real dt)
{
    vel_grad_[index_i] *= kernel_correction_(index_i);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
