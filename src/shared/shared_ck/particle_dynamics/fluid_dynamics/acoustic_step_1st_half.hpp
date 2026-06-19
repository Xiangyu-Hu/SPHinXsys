#ifndef ACOUSTIC_STEP_1ST_HALF_HPP
#define ACOUSTIC_STEP_1ST_HALF_HPP

#include "acoustic_step_1st_half.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
AcousticStep<BaseInteractionType>::AcousticStep(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
      dv_p_(this->particles_->template registerStateVariable<Real>("Pressure")),
      dv_compression_(this->particles_->template registerStateVariable<Real>("Compression", 1.0)),
      dv_compression_rate_(this->particles_->template registerStateVariable<Real>("CompressionRate")),
      dv_vel_(this->particles_->template registerStateVariable<Vecd>("Velocity")),
      dv_dpos_(this->particles_->template getVariableByName<Vecd>("Displacement")),
      dv_force_(this->particles_->template registerStateVariable<Vecd>("Force")),
      dv_force_prior_(this->particles_->template registerStateVariable<Vecd>("ForcePrior"))
{
    //----------------------------------------------------------------------
    //		add evolving variables
    //----------------------------------------------------------------------
    this->particles_->template addEvolvingVariable<Vecd>("Velocity");
    this->particles_->template addEvolvingVariable<Real>("Mass");
    this->particles_->template addEvolvingVariable<Vecd>("ForcePrior");
    this->particles_->template addEvolvingVariable<Real>("Compression");
    this->particles_->template addEvolvingVariable<Real>("CompressionRate");
    //----------------------------------------------------------------------
    //		add output particle data
    //----------------------------------------------------------------------
    this->particles_->template addVariableToWrite<Vecd>("Velocity");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class DynamicsIdentifier>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(DynamicsIdentifier &identifier)
    : AcousticStep<Interaction<Inner<Parameters...>>>(identifier),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_->getMatterMaterial())),
      riemann_solver_(this->fluid_, this->fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : eos_(ex_policy, encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedDataView(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataView(ex_policy)),
      compression_(encloser.dv_compression_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataView(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    compression_[index_i] += 0.5 * dt * compression_rate_[index_i];
    rho_[index_i] = compression_[index_i] * eos_.getReferenceDensity(index_i);
    p_[index_i] = eos_.PressureFromDensity(index_i, rho_[index_i]);
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_(ex_policy, encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataView(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    compression_rate_[index_i] = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        force_[index_i] -=
            riemann_.AverageP(
                index_i, index_j,
                static_cast<CorrectionDataType>(correction_(index_j) * p_[index_i]),
                static_cast<CorrectionDataType>(correction_(index_i) * p_[index_j])) *
            2.0 * dW_ijV_j * Vol_[index_i] * e_ij;
        compression_rate_[index_i] +=
            riemann_.DissipativeUJump(index_i, index_j, p_[index_i] - p_[index_j]) * dW_ijV_j;
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mass_(encloser.dv_mass_->DelegatedDataView(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataView(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataView(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class DynamicsIdentifier>
AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier), Interaction<Wall>(identifier),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_->getMatterMaterial())),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_(ex_policy, encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataView(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedDataView(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataView(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataView(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedDataView(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedDataView(ex_policy)),
      wall_acc_ave_(encloser.dv_wall_acc_ave_[contact_index]->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        Real face_wall_external_acceleration =
            (force_prior_[index_i] / mass_[index_i] - wall_acc_ave_[index_j]).dot(-e_ij);
        Real p_j_in_wall =
            p_[index_i] + rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
        force_[index_i] -=
            (p_[index_i] + p_j_in_wall) * correction_(index_i) * dW_ijV_j * Vol_[index_i] * e_ij;
        compression_rate_[index_i] +=
            riemann_.DissipativeUJump(index_i, index_j, p_[index_i] - p_j_in_wall) * dW_ijV_j;
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class DynamicsIdentifier>
AcousticStep1stHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier), kernel_correction_(this->particles_)
{
    SourceFluidType &source_fluid =
        DynamicCast<SourceFluidType>(this, this->sph_body_->getMatterMaterial());
    for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
    {
        contact_kernel_corrections_.push_back(KernelCorrectionType(this->contact_particles_[k]));
        TargetFluidType &target_fluid =
            DynamicCast<TargetFluidType>(this, this->contact_bodies_[k]->getMatterMaterial());
        riemann_solvers_.push_back(RiemannSolverType(source_fluid, target_fluid));
        dv_contact_p_.push_back(
            this->contact_particles_[k]->template registerStateVariable<Real>("Pressure"));
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      contact_correction_(ex_policy, encloser.contact_kernel_corrections_[contact_index]),
      riemann_(ex_policy, encloser.riemann_solvers_[contact_index]),
      Vol_(encloser.dv_Vol_->DelegatedDataView(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataView(ex_policy)),
      compression_rate_(encloser.dv_compression_rate_->DelegatedDataView(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataView(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedDataView(ex_policy)),
      contact_p_(encloser.dv_contact_p_[contact_index]->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        force_[index_i] -=
            riemann_.AverageP(
                index_i, index_j,
                static_cast<CorrectionDataType>(contact_correction_(index_j) * p_[index_i]),
                static_cast<CorrectionDataType>(correction_(index_i) * contact_p_[index_j])) *
            2.0 * dW_ijV_j * Vol_[index_i] * e_ij;
        compression_rate_[index_i] +=
            riemann_.DissipativeUJump(index_i, index_j, p_[index_i] - contact_p_[index_j]) * dW_ijV_j;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_1ST_HALF_HPP
