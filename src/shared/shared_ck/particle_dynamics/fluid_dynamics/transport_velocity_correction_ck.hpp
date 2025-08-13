#ifndef TRANSPORT_VELOCITY_CORRECTION_CK_HPP
#define TRANSPORT_VELOCITY_CORRECTION_CK_HPP

#include "transport_velocity_correction_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class KernelCorrectionType, class LimiterType, class ParticleScopeType, typename... Parameters>
TransportVelocityCorrectionCK<
    Inner<WithUpdate, KernelCorrectionType, LimiterType, ParticleScopeType, Parameters...>>::
    TransportVelocityCorrectionCK(Inner<Parameters...> &inner_relation, Real coefficient)
    : BaseInteraction(inner_relation),
      h_ref_(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
      correction_scaling_(coefficient * h_ref_ * h_ref_),
      h_ratio_(inner_relation.getNeighborMethod()), limiter_(h_ref_ * h_ref_),
      within_scope_method_(this->particles_),
      dv_dpos_(this->particles_->template getVariableByName<Vecd>("Displacement"))
{
    static_assert(std::is_base_of<WithinScope, ParticleScopeType>::value,
                  "WithinScope is not the base of ParticleScope!");
    static_assert(std::is_base_of<Limiter, LimiterType>::value,
                  "Limiter is not the base of LimiterType!");
}
//=================================================================================================//
template <class KernelCorrectionType, class LimiterType, class ParticleScopeType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
TransportVelocityCorrectionCK<
    Inner<WithUpdate, KernelCorrectionType, LimiterType, ParticleScopeType, Parameters...>>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : correction_scaling_(encloser.correction_scaling_),
      h_ratio_(encloser.h_ratio_), limiter_(encloser.limiter_),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)),
      kernel_gradient_integral_(encloser.dv_kernel_gradient_integral_->DelegatedData(ex_policy)),
      within_scope_(ex_policy, encloser.within_scope_method_, *this) {}
//=================================================================================================//
template <class KernelCorrectionType, class LimiterType, class ParticleScopeType, typename... Parameters>
void TransportVelocityCorrectionCK<
    Inner<WithUpdate, KernelCorrectionType, LimiterType, ParticleScopeType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        Real inv_h_ratio = 1.0 / h_ratio_(index_i);
        Vecd residual = this->kernel_gradient_integral_[index_i];
        Real squared_norm = residual.squaredNorm();
        dpos_[index_i] += correction_scaling_ * limiter_(squared_norm) *
                          this->kernel_gradient_integral_[index_i] * inv_h_ratio * inv_h_ratio;
    }
}
} // end namespace fluid_dynamics
} // end namespace SPH
#endif // TRANSPORT_VELOCITY_CORRECTION_CK_HPP
