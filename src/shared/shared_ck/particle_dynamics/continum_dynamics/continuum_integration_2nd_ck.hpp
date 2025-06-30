#ifndef CONTINUUM_INTERATION_2ND_CK_HPP
#define CONTINUUM_INTERATION_2ND_CK_HPP

#include "acoustic_step_2nd_half.h"

namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    PlasticAcousticStep2ndHalf(Inner<Parameters...> &inner_relation)
    : PlasticAcousticStep<Interaction<Inner<Parameters...>>>(inner_relation),
      correction_(this->particles_), riemann_solver_(this->plastic_continuum_, this->plastic_continuum_, 20.0 * (Real)Dimensions)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      velocity_gradient_(encloser.dv_velocity_gradient_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    Matd velocity_gradient = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = correction_(index_i) * this->e_ij(index_i, index_j);
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];

        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
        velocity_gradient -= (vel_[index_i] - vel_[index_j]) * dW_ijV_j * e_ij.transpose();
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] = p_dissipation * Vol_[index_i];
    velocity_gradient_[index_i] = velocity_gradient;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : plastic_kernel_(encloser.plastic_continuum_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      velocity_gradient_(encloser.dv_velocity_gradient_->DelegatedData(ex_policy)),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)),
      strain_tensor_3D_(encloser.dv_strain_tensor_3D_->DelegatedData(ex_policy)),
      stress_rate_3D_(encloser.dv_stress_rate_3D_->DelegatedData(ex_policy)),
      strain_rate_3D_(encloser.dv_strain_rate_3D_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    Mat3d velocity_gradient = upgradeToMat3d(velocity_gradient_[index_i]);
    Mat3d stress_tensor_rate_3D_ = plastic_kernel_.ConstitutiveRelation(velocity_gradient, stress_tensor_3D_[index_i]);
    stress_rate_3D_[index_i] += stress_tensor_rate_3D_; // stress diffusion is on
    stress_tensor_3D_[index_i] += stress_rate_3D_[index_i] * dt;
    /*return mapping*/
    stress_tensor_3D_[index_i] = plastic_kernel_.ReturnMapping(stress_tensor_3D_[index_i]);
    strain_rate_3D_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    strain_tensor_3D_[index_i] += strain_rate_3D_[index_i] * dt;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
PlasticAcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    PlasticAcousticStep2ndHalf(Contact<Parameters...> &wall_contact_relation)
    : BaseInteraction(wall_contact_relation), Interaction<Wall>(wall_contact_relation),
      correction_(this->particles_), riemann_solver_(this->plastic_continuum_, this->plastic_continuum_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
PlasticAcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedData(ex_policy)),
      wall_vel_ave_(encloser.dv_wall_vel_ave_[contact_index]->DelegatedData(ex_policy)),
      wall_n_(encloser.dv_wall_n_[contact_index]->DelegatedData(ex_policy)),
      velocity_gradient_(encloser.dv_velocity_gradient_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void PlasticAcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    Vecd vel_i = vel_[index_i];
    Matd velocity_gradient = Matd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * wall_Vol_[index_j];
        Vecd vel_in_wall = 2.0 * wall_vel_ave_[index_j] - vel_[index_i];
        density_change_rate += (vel_i - vel_in_wall).dot(e_ij) * dW_ijV_j;
        Real u_jump = 2.0 * (vel_i - wall_vel_ave_[index_j]).dot(wall_n_[index_j]);
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * wall_n_[index_j];
        velocity_gradient -= (vel_i - vel_in_wall) * dW_ijV_j * e_ij.transpose();
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation * Vol_[index_i];
    velocity_gradient_[index_i] += velocity_gradient;
}
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_INTERATION_2ND_CK_HPP