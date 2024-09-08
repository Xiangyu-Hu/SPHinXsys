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
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, this->sph_body_.getBaseMaterial())),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
      dv_p_(this->particles_->template registerStateVariableOnly<Real>("Pressure")),
      dv_drho_dt_(this->particles_->template registerStateVariableOnly<Real>("DensityChangeRate")),
      dv_vel_(this->particles_->template registerStateVariableOnly<Vecd>("Velocity")),
      dv_dpos_(this->particles_->template getVariableByName<Vecd>("Displacement")),
      dv_force_(this->particles_->template registerStateVariableOnly<Vecd>("Force")),
      dv_force_prior_(this->particles_->template registerStateVariableOnly<Vecd>("ForcePrior"))
{
    //----------------------------------------------------------------------
    //		add sortable particle data
    //----------------------------------------------------------------------
    this->particles_->template addVariableToSort<Vecd>("Position");
    this->particles_->template addVariableToSort<Vecd>("Velocity");
    this->particles_->template addVariableToSort<Real>("Mass");
    this->particles_->template addVariableToSort<Vecd>("ForcePrior");
    this->particles_->template addVariableToSort<Vecd>("Force");
    this->particles_->template addVariableToSort<Real>("DensityChangeRate");
    this->particles_->template addVariableToSort<Real>("Density");
    this->particles_->template addVariableToSort<Real>("Pressure");
    this->particles_->template addVariableToSort<Real>("VolumetricMeasure");
    //----------------------------------------------------------------------
    //		add restart output particle data
    //----------------------------------------------------------------------
    this->particles_->template addVariableToRestart<Vecd>("Position");
    this->particles_->template addVariableToRestart<Real>("VolumetricMeasure");
    this->particles_->template addVariableToRestart<Real>("Pressure");
    this->particles_->template addVariableToRestart<Real>("DensityChangeRate");
    this->particles_->template addVariableToRestart<Vecd>("Velocity");
    this->particles_->template addVariableToRestart<Vecd>("Force");
    this->particles_->template addVariableToRestart<Vecd>("ForcePrior");
    //----------------------------------------------------------------------
    //		add output particle data
    //----------------------------------------------------------------------
    this->particles_->template addVariableToWrite<Vecd>("Velocity");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(InnerRelation &inner_relation)
    : AcousticStep<Interaction<Inner<Parameters...>>>(inner_relation),
      correction_(this->particles_), riemann_solver_(this->fluid_, this->fluid_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrection is not the base of KernelCorrectionType!");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : eos_(encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    p_[index_i] = eos_.getPressure(rho_[index_i]);
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        force -= (p_[index_i] * correction_(index_j) + p_[index_j] * correction_(index_i)) * dW_ijV_j * e_ij;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
    }
    force_[index_i] += force * Vol_[index_i];
    drho_dt_[index_i] = rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(ContactRelation &wall_contact_relation)
    : AcousticStep<Interaction<Contact<Wall, Parameters...>>>(wall_contact_relation),
      correction_(this->particles_), riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(encloser.correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataField(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedDataField(ex_policy)),
      wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedDataField(ex_policy)),
      wall_acc_ave_(encloser.dv_wall_acc_ave_[contact_index]->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * wall_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        Real face_wall_external_acceleration = (force_prior_[index_i] / mass_[index_i] - wall_acc_ave_[index_j]).dot(-e_ij);
        Real p_in_wall = p_[index_i] + rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
        force -= (p_[index_i] + p_in_wall) * correction_(index_i) * dW_ijV_j * e_ij;
        rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_in_wall) * dW_ijV_j;
    }
    force_[index_i] += force * Vol_[index_i];
    drho_dt_[index_i] += rho_dissipation * rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_1ST_HALF_HPP
