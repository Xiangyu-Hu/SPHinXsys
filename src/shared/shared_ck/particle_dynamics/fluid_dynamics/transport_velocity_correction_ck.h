#ifndef TRANSPORT_VELOCITY_CORRECTION_CK_H
#define TRANSPORT_VELOCITY_CORRECTION_CK_H

#include "base_fluid_dynamics.h"
#include "particle_functors_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
template <class DynamicsIdentifier, class LimiterType, typename... ParticleScopes>
class TransportVelocityCorrectionCK : public BaseLocalDynamics<DynamicsIdentifier>
{
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopes...>::ComputingKernel;
    using Adaptation = typename DynamicsIdentifier::Adaptation;
    using SmoothingLengthRatio = typename Adaptation::SmoothingLengthRatioType;

  public:
    explicit TransportVelocityCorrectionCK(DynamicsIdentifier &identifier, Real coefficient = 0.2);
    virtual ~TransportVelocityCorrectionCK() {}

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real correction_scaling_;
        SmoothingLengthRatio h_ratio_;
        LimiterType limiter_;
        Vecd *dpos_;
        Vecd *kernel_gradient_integral_;
        ParticleScopeTypeKernel within_scope_;
    };

  protected:
    Real h_ref_;              ///< e.g. reference smoothing length
    Real correction_scaling_; ///< typically coefficient * h_ref^2
    LimiterType limiter_;     ///< e.g. a limiter on the final correction step
    ParticleScopeTypeCK<ParticleScopes...> within_scope_method_;
    DiscreteVariable<Vecd> *dv_dpos_, *dv_kernel_gradient_integral_;
    Adaptation &adaptation_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // TRANSPORT_VELOCITY_CORRECTION_CK_H
