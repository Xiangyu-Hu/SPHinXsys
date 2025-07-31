#ifndef EROSION_ACOUSTIC_STEP_2ND_HALF_HPP
#define EROSION_ACOUSTIC_STEP_2ND_HALF_HPP

#include "erosion_acoustic_step_2nd_half.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep2ndHalf<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(Contact<Parameters...> &soil_contact_relation)
    : BaseInteraction(soil_contact_relation), Interaction<Soil>(soil_contact_relation),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_.getBaseMaterial())),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      soil_Vol_(encloser.dv_soil_Vol_[contact_index]->DelegatedData(ex_policy)),
      soil_p_(encloser.dv_soil_p_[contact_index]->DelegatedData(ex_policy)),
      soil_n_(encloser.dv_soil_n_[contact_index]->DelegatedData(ex_policy)),
      soil_vel_(encloser.dv_soil_vel_[contact_index]->DelegatedData(ex_policy)){}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * soil_Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);
        Vecd vel_in_soil = 2.0 * soil_vel_[index_j]- vel_[index_i];
        density_change_rate += (vel_[index_i] - vel_in_soil).dot(corrected_e_ij) * dW_ijV_j;
        Real u_jump = 2.0 * (vel_[index_i] - soil_vel_[index_j]).dot(soil_n_[index_j]);
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * soil_n_[index_j]; 
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation * Vol_[index_i];
}
} // namespace fluid_dynamics
} // namespace SPH
#endif // EROSION_ACOUSTIC_STEP_2ND_HALF_HPP
