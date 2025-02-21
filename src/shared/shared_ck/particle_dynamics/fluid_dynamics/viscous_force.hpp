#pragma once

#include "viscous_force.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType,
          template <typename...> class RelationType, typename... Parameters>
template <class BaseRelationType>
ViscousForceCK<Base, ViscosityType, KernelCorrectionType, RelationType<Parameters...>>::
    ViscousForceCK(BaseRelationType &base_relation)
    : Interaction<RelationType<Parameters...>>(base_relation),
      ForcePriorCK(this->particles_, "ViscousForce"),
      viscosity_model_(DynamicCast<ViscosityType>(this, this->particles_->getBaseMaterial())),
      kernel_correction_(this->particles_),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      dv_viscous_force_(this->dv_current_force_),
      smoothing_length_sq_(pow(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength(), 2)) {}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType,
          template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
ViscousForceCK<Base, ViscosityType, KernelCorrectionType, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : Interaction<RelationType<Parameters...>>::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      viscosity_(ex_policy, encloser.viscosity_model_),
      correction_(ex_policy, encloser.kernel_correction_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      viscous_force_(encloser.dv_viscous_force_->DelegatedData(ex_policy)),
      smoothing_length_sq_(encloser.smoothing_length_sq_) {}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType, typename... Parameters>
void ViscousForceCK<Inner<WithUpdate, ViscosityType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);
        Vecd vel_derivative = (this->vel_[index_i] - this->vel_[index_j]) /
                              (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);

        force += vec_r_ij.dot((this->correction_(index_i) + this->correction_(index_j)) * e_ij) *
                 harmonic_average(this->viscosity_(index_i), this->viscosity_(index_j)) *
                 vel_derivative * this->dW_ij(index_i, index_j) * this->Vol_[index_j];
    }
    this->viscous_force_[index_i] = force * this->Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType, typename... Parameters>
void ViscousForceCK<Contact<Wall, ViscosityType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);
        Vecd vel_derivative = 2.0 * (this->vel_[index_i] - this->wall_vel_ave_[index_j]) /
                              (vec_r_ij.squaredNorm() + 0.01 * this->smoothing_length_sq_);

        force += 2.0 * vec_r_ij.dot(this->correction_(index_i) * e_ij) *
                 this->viscosity_(index_i) * vel_derivative *
                 this->dW_ij(index_i, index_j) * this->wall_Vol_[index_j];
    }
    this->viscous_force_[index_i] += force * this->Vol_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
