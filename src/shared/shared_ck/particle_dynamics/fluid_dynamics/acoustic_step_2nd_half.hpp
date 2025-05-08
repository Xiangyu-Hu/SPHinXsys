#ifndef ACOUSTIC_STEP_2ND_HALF_HPP
#define ACOUSTIC_STEP_2ND_HALF_HPP

#include "acoustic_step_2nd_half.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(Inner<Parameters...> &inner_relation)
    : AcousticStep<Interaction<Inner<Parameters...>>>(inner_relation),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_.getBaseMaterial())),
      riemann_solver_(this->fluid_, this->fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);

        Real u_jump = (vel_[index_i] - vel_[index_j]).dot(corrected_e_ij);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * corrected_e_ij;
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] = p_dissipation * Vol_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(Contact<Parameters...> &wall_contact_relation)
    : BaseInteraction(wall_contact_relation), Interaction<Wall>(wall_contact_relation),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_.getBaseMaterial())),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
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
      wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedData(ex_policy)),
      wall_vel_ave_(encloser.dv_wall_vel_ave_[contact_index]->DelegatedData(ex_policy)),
      wall_n_(encloser.dv_wall_n_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * wall_Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);

        Vecd vel_j_in_wall = 2.0 * wall_vel_ave_[index_j] - vel_[index_i];
        density_change_rate += (vel_[index_i] - vel_j_in_wall).dot(corrected_e_ij) * dW_ijV_j;
        Real u_jump = 2.0 * (vel_[index_i] - wall_vel_ave_[index_j]).dot(wall_n_[index_j]);
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * wall_n_[index_j];
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation * Vol_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep2ndHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep2ndHalf(Contact<Parameters...> &wall_contact_relation)
    : BaseInteraction(wall_contact_relation), kernel_correction_(this->particles_)
{
    SourceFluidType &source_fluid =
        DynamicCast<SourceFluidType>(this, this->sph_body_.getBaseMaterial());
    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        TargetFluidType &target_fluid =
            DynamicCast<TargetFluidType>(this, this->contact_bodies_[k]->getBaseMaterial());
        riemann_solvers_.push_back(RiemannSolverType(source_fluid, target_fluid));
        dv_contact_Vol_.push_back(
            this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
        dv_contact_vel_.push_back(
            this->contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep2ndHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_solver_(encloser.riemann_solvers_[contact_index]),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_vel_(encloser.dv_contact_vel_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep2ndHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        Vecd corrected_e_ij = correction_(index_i) * this->e_ij(index_i, index_j);

        Real u_jump = (vel_[index_i] - contact_vel_[index_j]).dot(corrected_e_ij);
        density_change_rate += u_jump * dW_ijV_j;
        p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * corrected_e_ij;
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    force_[index_i] += p_dissipation * Vol_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_2ND_HALF_HPP
