#ifndef EROSION_ACOUSTIC_STEP_1ST_HALF_HPP
#define EROSION_ACOUSTIC_STEP_1ST_HALF_HPP

#include "erosion_acoustic_step_1st_half.h"

namespace SPH
{
namespace fluid_dynamics
{
    //=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep1stHalf<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(Contact<Parameters...> &soil_contact_relation)
    : BaseInteraction(soil_contact_relation), Interaction<Soil>(soil_contact_relation),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_.getBaseMaterial())),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      soil_Vol_(encloser.dv_soil_Vol_[contact_index]->DelegatedData(ex_policy)),
      soil_p_(encloser.dv_soil_p_[contact_index]->DelegatedData(ex_policy)),
      soil_n_(encloser.dv_soil_n_[contact_index]->DelegatedData(ex_policy)),
      soil_vel_(encloser.dv_soil_vel_[contact_index]->DelegatedData(ex_policy)),
      soil_shear_stress_tensor_(encloser.dv_soil_shear_stress_tensor_[contact_index]->DelegatedData(ex_policy)){}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Contact<Soil, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    Matd soil_shear_stress_tensor_i = degradeToMatd(Mat3d::Zero());
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * soil_Vol_[index_j];
        Vecd nablaW_ijV_j = this->dW_ij(index_i, index_j) * soil_Vol_[index_j] * this->e_ij(index_i, index_j);
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();
        Real face_wall_external_acceleration = (force_prior_[index_i] / mass_[index_i]).dot(-e_ij);
        Real p_j_in_soil = p_[index_i] + rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
        force -= (p_[index_i] + p_j_in_soil) * correction_(index_i) * dW_ijV_j * e_ij;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_j_in_soil) * dW_ijV_j;

        Matd soil_shear_stress_tensor_j = degradeToMatd(soil_shear_stress_tensor_[index_j]);
        force += ((soil_shear_stress_tensor_i + soil_shear_stress_tensor_j)) * nablaW_ijV_j;
    }
    force_[index_i] += force * Vol_[index_i];
    drho_dt_[index_i] += rho_dissipation * rho_[index_i];
}
} // namespace fluid_dynamics
} // namespace SPH
#endif // EROSION_ACOUSTIC_STEP_1ST_HALF_HPP
