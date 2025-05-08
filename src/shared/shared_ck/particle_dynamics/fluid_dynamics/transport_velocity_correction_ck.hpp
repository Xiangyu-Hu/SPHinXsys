#ifndef TRANSPORT_VELOCITY_CORRECTION_CK_HPP
#define TRANSPORT_VELOCITY_CORRECTION_CK_HPP
#include "transport_velocity_correction_ck.h"

namespace SPH
{
namespace fluid_dynamics
{

//========================================================================================
//   1) Implementation of the Base Class
//========================================================================================
template <class BaseInteractionType>
template <class DynamicsIdentifier>
TransportVelocityCorrectionCKBase<BaseInteractionType>::
    TransportVelocityCorrectionCKBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_dpos_(this->particles_->template getVariableByName<Vecd>("Displacement")),
      dv_zero_gradient_residue_(
          this->particles_->template registerStateVariableOnly<Vecd>("ZeroGradientResidue"))
{
}

//========================================================================================
//   2) Partial Specialization:
//      <Inner<WithUpdate, KernelCorrectionType, ResolutionType, LimiterType, ParticleScope, ExtraParams...>
//========================================================================================
template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
TransportVelocityCorrectionCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>::
    TransportVelocityCorrectionCK(Inner<Parameters...> &inner_relation, Real coefficient)
    : TransportVelocityCorrectionCKBase<Interaction<Inner<Parameters...>>>(inner_relation),
      kernel_correction_(this->particles_),
      h_ref_(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
      correction_scaling_(coefficient * h_ref_ * h_ref_),
      h_ratio_(this->particles_),
      limiter_(h_ref_ * h_ref_),
      within_scope_method_(this->particles_)
{
    static_assert(std::is_base_of<WithinScope, ParticleScopeType>::value,
                  "WithinScope is not the base of ParticleScope!");
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrectionType must derive from KernelCorrection!");
    static_assert(std::is_base_of<Limiter, LimiterType>::value,
                  "Limiter is not the base of LimiterType!");
}

//------------------------------------------------------------------------------
// InteractKernel Implementation
//------------------------------------------------------------------------------
template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
TransportVelocityCorrectionCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)),
      zero_gradient_residue_(encloser.dv_zero_gradient_residue_->DelegatedData(ex_policy)),
      within_scope_(ex_policy, encloser.within_scope_method_, *this)
{
}

template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
void TransportVelocityCorrectionCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    {
        Vecd inconsistency = Vecd::Zero();
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            UnsignedInt index_j = this->neighbor_index_[n];
            const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
            const Vecd e_ij = this->e_ij(index_i, index_j);

            // acceleration for transport velocity
            inconsistency -= (this->correction_(index_i) + this->correction_(index_j)) *
                             dW_ijV_j * e_ij;
        }
        this->zero_gradient_residue_[index_i] = inconsistency;
    }
}

//------------------------------------------------------------------------------
// UpdateKernel Implementation
//------------------------------------------------------------------------------
template <class UpdatePolicy,
          class KernelCorrectionType,
          class ResolutionType,
          class LimiterType,
          class ParticleScopeType,
          typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
TransportVelocityCorrectionCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : correction_(ex_policy, encloser.kernel_correction_),
      correction_scaling_(encloser.correction_scaling_),
      h_ratio_(encloser.h_ratio_),
      limiter_(encloser.limiter_),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)),
      zero_gradient_residue_(encloser.dv_zero_gradient_residue_->DelegatedData(ex_policy)), within_scope_(ex_policy, encloser.within_scope_method_, *this)
{
}

template <class UpdatePolicy,
          class KernelCorrectionType,
          class ResolutionType,
          class LimiterType,
          class ParticleScopeType,
          typename... ExtraParams>
void TransportVelocityCorrectionCK<
    Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, ExtraParams...>>::UpdateKernel::
    update(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        Real inv_h_ratio = 1.0 / h_ratio_(index_i);
        Vecd residue = this->zero_gradient_residue_[index_i];
        Real squared_norm = residue.squaredNorm();
        dpos_[index_i] += correction_scaling_ * limiter_(squared_norm) *
                          this->zero_gradient_residue_[index_i] * inv_h_ratio * inv_h_ratio;
    }
}

//========================================================================================
//   Partial Specialization:
//      <Contact<WithUpdate, KernelCorrectionType, ParticleScopeType, ExtraParams...>>
//========================================================================================

template <class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
TransportVelocityCorrectionCK<Contact<Wall, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>::
    TransportVelocityCorrectionCK(Contact<Parameters...> &contact_relation)
    : BaseInteraction(contact_relation), Interaction<Wall>(contact_relation),
      kernel_correction_(this->particles_)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        dv_contact_wall_Vol_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//------------------------------------------------------------------------------
// InteractKernel Implementation
//------------------------------------------------------------------------------
template <class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
TransportVelocityCorrectionCK<Contact<Wall, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>::
    InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      zero_gradient_residue_(encloser.dv_zero_gradient_residue_->DelegatedData(ex_policy)),
      contact_wall_Vol_(encloser.dv_contact_wall_Vol_[contact_index]->DelegatedData(ex_policy))
{
}

template <class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
void TransportVelocityCorrectionCK<Contact<Wall, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd inconsistency = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_wall_Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);

        // acceleration for transport velocity
        inconsistency -= 2.0 * this->correction_(index_i) *
                         dW_ijV_j * e_ij;
    }
    this->zero_gradient_residue_[index_i] += inconsistency;
}

} // end namespace fluid_dynamics
} // end namespace SPH
#endif // TRANSPORT_VELOCITY_CORRECTION_CK_HPP
