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
      dv_drho_dt_(this->particles_->template registerStateVariable<Real>("DensityChangeRate")),
      dv_vel_(this->particles_->template registerStateVariable<Vecd>("Velocity")),
      dv_dpos_(this->particles_->template getVariableByName<Vecd>("Displacement")),
      dv_force_(this->particles_->template registerStateVariable<Vecd>("Force")),
      dv_force_prior_(this->particles_->template registerStateVariable<Vecd>("ForcePrior")),
      dv_p_background_(this->particles_->template registerStateVariable<Real>("BackgroundPressure")),
      dv_perturbation_momentum_(this->particles_->template registerStateVariable<Real>("PerturbationMomentum")),
      dv_acoustic_dens_mode_(this->particles_->template registerStateVariable<Real>("AcousticDensMode")),
      dv_reference_body_acceleration_(
          this->particles_->template registerStateVariable<Vecd>("ReferenceBodyAcceleration"))
{
    //----------------------------------------------------------------------
    //		add evolving variables
    //----------------------------------------------------------------------
    this->particles_->template addEvolvingVariable<Vecd>("Velocity");
    this->particles_->template addEvolvingVariable<Real>("Mass");
    this->particles_->template addEvolvingVariable<Vecd>("ForcePrior");
    this->particles_->template addEvolvingVariable<Real>("DensityChangeRate");
    this->particles_->template addEvolvingVariable<Real>("Density");
    this->particles_->template addEvolvingVariable<Real>("Pressure");
    this->particles_->template addEvolvingVariable<Real>("BackgroundPressure");
    this->particles_->template addEvolvingVariable<Real>("PerturbationMomentum");
    this->particles_->template addEvolvingVariable<Real>("AcousticDensMode");
    this->particles_->template addEvolvingVariable<Vecd>("ReferenceBodyAcceleration");
    //----------------------------------------------------------------------
    //		add output particle data
    //----------------------------------------------------------------------
    this->particles_->template addVariableToWrite<Vecd>("Velocity");
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(Inner<Parameters...> &inner_relation)
    : AcousticStep<Interaction<Inner<Parameters...>>>(inner_relation),
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
    : eos_(encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      p_background_(encloser.dv_p_background_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InitializeKernel::initialize(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
    // Doc87 R3-A: p = p' + p_h with p' = EOS(ρ). Default p_h=0 → legacy behavior.
    p_[index_i] = eos_.getPressure(rho_[index_i]) + p_background_[index_i];
    dpos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      p_background_(encloser.dv_p_background_->DelegatedData(ex_policy)),
      perturbation_momentum_(encloser.dv_perturbation_momentum_->DelegatedData(ex_policy)),
      acoustic_dens_mode_(encloser.dv_acoustic_dens_mode_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    // Doc105: dens mode 0 = legacy total-p; mode 1/2 = p'; mode 3 = omit p dens jump.
    // Momentum still via PerturbationMomentum.
    const Real dens_mode = acoustic_dens_mode_[index_i];
    const Real alpha = perturbation_momentum_[index_i];
    const Real p_i_prime = p_[index_i] - p_background_[index_i];
    const Real p_i_mom = p_[index_i] - alpha * p_background_[index_i];
    const bool dens_use_prime = dens_mode > Real(0.5) && dens_mode <= Real(2.5);
    const bool dens_omit_p_jump = dens_mode > Real(2.5);
    const Real p_i_dens = dens_use_prime ? p_i_prime : p_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        const Real p_j_prime = p_[index_j] - p_background_[index_j];
        const Real p_j_mom = p_[index_j] - alpha * p_background_[index_j];
        const Real p_j_dens = dens_use_prime ? p_j_prime : p_[index_j];

        force -= riemann_solver_.AverageP(
                     index_i, index_j,
                     static_cast<CorrectionDataType>(correction_(index_j) * p_i_mom),
                     static_cast<CorrectionDataType>(correction_(index_i) * p_j_mom)) *
                 2.0 * dW_ijV_j * e_ij;
        if (!dens_omit_p_jump)
            rho_dissipation += riemann_solver_.DissipativeUJump(index_i, index_j, p_i_dens - p_j_dens) * dW_ijV_j;
    }
    force_[index_i] += force * Vol_[index_i];
    drho_dt_[index_i] = rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)) {}
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
    AcousticStep1stHalf(Contact<Parameters...> &wall_contact_relation)
    : BaseInteraction(wall_contact_relation), Interaction<Wall>(wall_contact_relation),
      kernel_correction_(this->particles_),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_->getMatterMaterial())),
      riemann_solver_(this->fluid_, this->fluid_)
{
    // Doc201 B2 / Doc202: share wall tags with 2nd-half (idempotent register).
    const StdVec<BaseParticles *> contact_particles = wall_contact_relation.getContactParticles();
    for (size_t k = 0; k != contact_particles.size(); ++k)
    {
        BaseParticles *wp = contact_particles[k];
        dv_wall_free_slip_.push_back(wp->template registerDiscreteVariable<Real>(
            "AcousticWallFreeSlip", wp->ParticlesBound(), Real(0)));
        dv_wall_contact_omit_.push_back(wp->template registerDiscreteVariable<Real>(
            "AcousticWallContactOmit", wp->ParticlesBound(), Real(0)));
        dv_wall_contact_omit_1st_.push_back(wp->template registerDiscreteVariable<Real>(
            "AcousticWallContactOmit1st", wp->ParticlesBound(), Real(0)));
        // Keep 2nd-half fields registered so host tags can always write.
        wp->template registerDiscreteVariable<Real>(
            "AcousticWallContactOmit2nd", wp->ParticlesBound(), Real(0));
        wp->template registerDiscreteVariable<Real>(
            "AcousticWallContSilent", wp->ParticlesBound(), Real(0));
        wp->template registerDiscreteVariable<Real>(
            "AcousticWallOmitPJump", wp->ParticlesBound(), Real(0));
    }
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::InteractKernel(
        const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      riemann_solver_(encloser.riemann_solver_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      p_background_(encloser.dv_p_background_->DelegatedData(ex_policy)),
      perturbation_momentum_(encloser.dv_perturbation_momentum_->DelegatedData(ex_policy)),
      acoustic_dens_mode_(encloser.dv_acoustic_dens_mode_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      reference_body_acceleration_(encloser.dv_reference_body_acceleration_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      wall_acc_ave_(encloser.dv_wall_acc_ave_[contact_index]->DelegatedData(ex_policy)),
      wall_n_(encloser.dv_wall_n_[contact_index]->DelegatedData(ex_policy)),
      wall_free_slip_(encloser.dv_wall_free_slip_[contact_index]->DelegatedData(ex_policy)),
      wall_contact_omit_(encloser.dv_wall_contact_omit_[contact_index]->DelegatedData(ex_policy)),
      wall_contact_omit_1st_(encloser.dv_wall_contact_omit_1st_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    const Real dens_mode = acoustic_dens_mode_[index_i];
    const Real alpha = perturbation_momentum_[index_i];
    const Real p_i_prime = p_[index_i] - p_background_[index_i];
    const Real p_i_mom = p_[index_i] - alpha * p_background_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        // Doc202: full omit (H1) or 1st-half-only omit (H1c-a).
        if (wall_contact_omit_[index_j] > Real(0.5) || wall_contact_omit_1st_[index_j] > Real(0.5))
            continue;
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);
        Real r_ij = this->vec_r_ij(index_i, index_j).norm();

        // Actual external acceleration used for Adami ghost (may include gravity, buoyancy, ...).
        Real face_wall_external_acceleration =
            (force_prior_[index_i] / mass_[index_i] - wall_acc_ave_[index_j]).dot(-e_ij);
        const Real C_ext = rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);

        // Momentum: Adami ghost on the momentum pressure variable (total p or p').
        // Mode 2 (Doc186) only: omit C_ext when body-force hydro is already in p'.
        // Mode 3 (Doc201) must keep C_ext (same mom as mode 1); only dens jump is omitted below.
        const bool mom_omit_c_ext = dens_mode > Real(1.5) && dens_mode <= Real(2.5);
        Real p_j_mom = p_i_mom + (mom_omit_c_ext ? Real(0) : C_ext);
        Vecd f_pair = (p_i_mom + p_j_mom) * correction_(index_i) * dW_ijV_j * e_ij;
        // Doc201 B2 free-slip lid: pressure may not exert tangential traction (project onto n).
        // Legacy no-slip walls keep full e_ij force (Adami).
        if (wall_free_slip_[index_j] > Real(0.5))
        {
            const Vecd n_face = SGN(e_ij.dot(wall_n_[index_j])) * wall_n_[index_j];
            f_pair = f_pair.dot(n_face) * n_face;
        }
        force -= f_pair;

        /**
         * Doc105 wall dens:
         * mode 0 (legacy): DissipativeUJump(p_i_mom − p_j_mom) = −C_ext  (non-zero when C≠0).
         * mode 1 (p' dens): jump = p_i' − (p_i+C_ext − (p_h_i+C_ref)) = C_ext − C_ref + (p' terms cancel at rest).
         *   C_ref from ReferenceBodyAcceleration (not ForcePrior).
         *   P2b (g not in ForcePrior): set a_ref=0 → C_ref=0; do not set a_ref=g or dens sees +C_g.
         * mode 2 (Doc186): dens same as mode 1 (still uses C_ext/C_ref); momentum omits C_ext above.
         * mode 3 (Doc201): omit pressure dens jump (continuity from 2nd-half u only).
         * Doc201 B2 H-dens-fs: AcousticWallFreeSlip wall particles also omit DissipativeUJump
         *   (local mode-3 on free-slip faces only; side/bottom noslip dens unchanged).
         */
        const bool omit_wall_dens_jump =
            (dens_mode > Real(2.5)) || (wall_free_slip_[index_j] > Real(0.5));
        if (!omit_wall_dens_jump)
        {
            Real p_jump_for_dens = Real(0);
            if (dens_mode > Real(0.5))
            {
                const Real face_wall_reference_acceleration =
                    (reference_body_acceleration_[index_i] - wall_acc_ave_[index_j]).dot(-e_ij);
                const Real C_ref = rho_[index_i] * r_ij * SMAX(Real(0), face_wall_reference_acceleration);
                const Real p_j_total = p_[index_i] + C_ext;
                const Real p_h_j = p_background_[index_i] + C_ref;
                const Real p_j_prime = p_j_total - p_h_j;
                p_jump_for_dens = p_i_prime - p_j_prime;
            }
            else
            {
                p_jump_for_dens = p_i_mom - p_j_mom; // legacy Adami dens on momentum pressure
            }
            rho_dissipation += riemann_solver_.DissipativeUJump(index_i, index_j, p_jump_for_dens) * dW_ijV_j;
        }
    }
    force_[index_i] += force * Vol_[index_i];
    drho_dt_[index_i] += rho_dissipation * rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
AcousticStep1stHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    AcousticStep1stHalf(Contact<Parameters...> &contact_relation)
    : BaseInteraction(contact_relation), kernel_correction_(this->particles_)
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
      riemann_solver_(encloser.riemann_solvers_[contact_index]),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      drho_dt_(encloser.dv_drho_dt_->DelegatedData(ex_policy)),
      force_(encloser.dv_force_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_p_(encloser.dv_contact_p_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
void AcousticStep1stHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Real rho_dissipation(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        Vecd e_ij = this->e_ij(index_i, index_j);

        force -= riemann_solver_.AverageP(
                     index_i, index_j,
                     static_cast<CorrectionDataType>(contact_correction_(index_j) * p_[index_i]),
                     static_cast<CorrectionDataType>(correction_(index_i) * contact_p_[index_j])) *
                 2.0 * dW_ijV_j * e_ij;
        rho_dissipation += riemann_solver_.DissipativeUJump(index_i, index_j, p_[index_i] - contact_p_[index_j]) * dW_ijV_j;
    }
    force_[index_i] += force * Vol_[index_i];
    drho_dt_[index_i] += rho_dissipation * rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_1ST_HALF_HPP
