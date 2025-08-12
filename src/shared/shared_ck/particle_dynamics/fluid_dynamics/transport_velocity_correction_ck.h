#ifndef TRANSPORT_VELOCITY_CORRECTION_CK_H
#define TRANSPORT_VELOCITY_CORRECTION_CK_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"
#include "kernel_gradient_integral.hpp"
#include "particle_functors_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename...>
class TransportVelocityCorrectionCK;

template <class KernelCorrectionType, class LimiterType, class ParticleScopeType, typename... Parameters>
class TransportVelocityCorrectionCK<
    Inner<WithUpdate, KernelCorrectionType, LimiterType, ParticleScopeType, Parameters...>>
    : public KernelGradientIntegral<Inner<KernelCorrectionType, Parameters...>>
{
    using BaseInteraction = KernelGradientIntegral<Inner<KernelCorrectionType, Parameters...>>;
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;
    using SmoothingRatio = typename Inner<Parameters...>::NeighborMethodType::SmoothingRatio;

  public:
    explicit TransportVelocityCorrectionCK(Inner<Parameters...> &inner_relation, Real coefficient = 0.2);
    virtual ~TransportVelocityCorrectionCK() {}

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real correction_scaling_;
        SmoothingRatio h_ratio_;
        LimiterType limiter_;
        Vecd *dpos_;
        Vecd *kernel_gradient_integral_;
        ParticleScopeTypeKernel within_scope_;
    };

  protected:
    Real h_ref_;              ///< e.g. reference smoothing length
    Real correction_scaling_; ///< typically coefficient * h_ref^2
    SmoothingRatio h_ratio_;  ///< e.g. for adaptive resolution
    LimiterType limiter_;     ///< e.g. a limiter on the final correction step
    ParticleScopeTypeCK<ParticleScopeType> within_scope_method_;
    DiscreteVariable<Vecd> *dv_dpos_;
};

template <class KernelCorrectionType, typename... Parameters>
class TransportVelocityCorrectionCK<Contact<Boundary, KernelCorrectionType, Parameters...>>
    : public KernelGradientIntegral<Contact<Boundary, KernelCorrectionType, Parameters...>>
{
    using BaseInteraction = KernelGradientIntegral<Contact<Boundary, KernelCorrectionType, Parameters...>>;

  public:
    explicit TransportVelocityCorrectionCK(Contact<Parameters...> &contact_relation)
        : BaseInteraction(contact_relation) {}
    virtual ~TransportVelocityCorrectionCK() {}
};
//--------------------------------------------------------------------------------------
// Alias Definitions for Specific Configurations
//--------------------------------------------------------------------------------------
using TransportVelocityCorrectionComplexBulkParticlesCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, NoKernelCorrectionCK, NoLimiter, BulkParticles>,
        Contact<Boundary, NoKernelCorrectionCK>>;

using TransportVelocityLimitedCorrectionCorrectedComplexBulkParticlesCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, LinearCorrectionCK, TruncatedLinear, BulkParticles>,
        Contact<Boundary, LinearCorrectionCK>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // TRANSPORT_VELOCITY_CORRECTION_CK_H
