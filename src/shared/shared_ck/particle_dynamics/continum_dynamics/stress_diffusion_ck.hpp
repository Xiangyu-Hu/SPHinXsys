#ifndef STRESS_DIFFUSION_CK_HPP
#define STRESS_DIFFUSION_CK_HPP

#include "stress_diffusion_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <typename... Parameters>
StressDiffusionCK<Inner<Parameters...>>::StressDiffusionCK(Inner<Parameters...> &inner_relation)
    : PlasticAcousticStep<Interaction<Inner<Parameters...>>>(inner_relation),
      dv_phi_(DynamicCast<PlasticContinuum>(this, this->plastic_continuum_).getFrictionAngle()),
      dv_smoothing_length_(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
      dv_sound_speed_(this->plastic_continuum_.ReferenceSoundSpeed()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")){};
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
StressDiffusionCK<Inner<Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      plastic_kernel_(encloser.plastic_continuum_),
      zeta_(encloser.dv_zeta_), phi_(encloser.dv_phi_),
      smoothing_length_(encloser.dv_smoothing_length_), sound_speed_(encloser.dv_sound_speed_),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)), Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
      stress_rate_3D_(encloser.dv_stress_rate_3D_->DelegatedData(ex_policy)){};
//=================================================================================================//
template <typename... Parameters>
void StressDiffusionCK<Inner<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd acc_prior_i = force_prior_[index_i] / this->mass_[index_i];
    Real gravity = abs(acc_prior_i(1, 0));
    Real density = plastic_kernel_.getDensity();
    Mat3d diffusion_stress_rate = Mat3d::Zero();
    Mat3d diffusion_stress = Mat3d::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Real y_ij = pos_[index_i](1, 0) - pos_[index_j](1, 0);
        diffusion_stress = stress_tensor_3D_[index_i] - stress_tensor_3D_[index_j];
        diffusion_stress(0, 0) -= (1 - math::sin(phi_)) * density * gravity * y_ij;
        diffusion_stress(1, 1) -= density * gravity * y_ij;
        diffusion_stress(2, 2) -= (1 - math::sin(phi_)) * density * gravity * y_ij;
        diffusion_stress_rate += 2 * zeta_ * smoothing_length_ * sound_speed_ *
                                 diffusion_stress * r_ij * dW_ijV_j / (r_ij * r_ij + 0.01 * smoothing_length_);
    }
    stress_rate_3D_[index_i] = diffusion_stress_rate;
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // STRESS_DIFFUSION_CK_HPP
