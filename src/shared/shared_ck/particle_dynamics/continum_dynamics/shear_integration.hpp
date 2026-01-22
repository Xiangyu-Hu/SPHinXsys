#ifndef SHEAR_INTEGRATION_HPP
#define SHEAR_INTEGRATION_HPP

#include "shear_integration.h"

#include "base_particles.hpp"

namespace SPH
{
namespace continuum_dynamics
{
//====================================================================================//
template <class MaterialType, typename... Parameters>
ShearIntegration<Inner<OneLevel, MaterialType, Parameters...>>::
    ShearIntegration(Inner<Parameters...> &inner_relation, Real xi)
    : BaseInteraction(inner_relation), ForcePriorCK(this->particles_, "ShearForce"),
      materal_(DynamicCast<MaterialType>(this, this->sph_body_->getBaseMaterial())),
      xi_(xi), dv_shear_force_(this->getCurrentForce()),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      dv_vel_gradient_(this->particles_->template getVariableByName<Matd>("VelocityGradient")),
      dv_strain_tensor_(this->particles_->template registerStateVariable<Matd>("StrainTensor")),
      dv_shear_stress_(this->particles_->template registerStateVariable<Matd>("ShearStress")),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_scale_penalty_force_(this->particles_->template registerStateVariable<Real>("ScalePenaltyForce"))
{
    this->particles_->template addEvolvingVariable<Matd>(dv_strain_tensor_);
    this->particles_->template addEvolvingVariable<Matd>(dv_shear_stress_);
}
//====================================================================================//
template <class MaterialType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ShearIntegration<Inner<OneLevel, MaterialType, Parameters...>>::InitializeKernel::
    InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : constitute_(ex_policy, encloser.materal_), xi_(encloser.xi_),
      vel_gradient_(encloser.dv_vel_gradient_->DelegatedData(ex_policy)),
      strain_tensor_(encloser.dv_strain_tensor_->DelegatedData(ex_policy)),
      shear_stress_(encloser.dv_shear_stress_->DelegatedData(ex_policy)),
      scale_penalty_force_(encloser.dv_scale_penalty_force_->DelegatedData(ex_policy)) {}
//====================================================================================//
template <class MaterialType, typename... Parameters>
void ShearIntegration<Inner<OneLevel, MaterialType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    Matd strain_rate = 0.5 * (vel_gradient_[index_i] + vel_gradient_[index_i].transpose());
    strain_tensor_[index_i] += strain_rate * dt;

    Matd shear_stress_rate =
        constitute_.ShearStressRate(index_i, vel_gradient_[index_i], shear_stress_[index_i]);
    Matd shear_stress_try = shear_stress_[index_i] + shear_stress_rate * dt;
    scale_penalty_force_[index_i] = xi_ * constitute_.ScalePenaltyForce(index_i, shear_stress_try);
    shear_stress_[index_i] = constitute_.updateShearStress(index_i, shear_stress_try);
}
//====================================================================================//
template <class MaterialType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ShearIntegration<Inner<OneLevel, MaterialType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      G_(encloser.materal_.ShearModulus()),
      shear_force_(encloser.dv_shear_force_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      vel_gradient_(encloser.dv_vel_gradient_->DelegatedData(ex_policy)),
      shear_stress_(encloser.dv_shear_stress_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      scale_penalty_force_(encloser.dv_scale_penalty_force_->DelegatedData(ex_policy)) {}
//====================================================================================//
template <class MaterialType, typename... Parameters>
void ShearIntegration<Inner<OneLevel, MaterialType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd sum = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Vecd vec_r_ij = this->vec_r_ij(index_i, index_j);

        sum += (shear_stress_[index_i] + shear_stress_[index_j]) * dW_ijV_j * e_ij;
        Vecd v_ij_correction = vel_[index_i] - vel_[index_j] -
                               0.5 * (vel_gradient_[index_i] + vel_gradient_[index_j]) * vec_r_ij;
        Real penalty_scale = 0.5 * (scale_penalty_force_[index_i] + scale_penalty_force_[index_j]);
        sum += penalty_scale * G_ * v_ij_correction.dot(e_ij) * vec_r_ij *
               dW_ijV_j * dt / vec_r_ij.squaredNorm();
    }
    shear_force_[index_i] = sum * Vol_[index_i];
}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
#endif // SHEAR_INTEGRATION_HPP