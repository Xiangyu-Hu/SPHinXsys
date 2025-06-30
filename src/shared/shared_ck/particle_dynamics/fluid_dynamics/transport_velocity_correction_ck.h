#ifndef TRANSPORT_VELOCITY_CORRECTION_CK_H
#define TRANSPORT_VELOCITY_CORRECTION_CK_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"
#include "particle_functors_ck.h"

namespace SPH
{
namespace fluid_dynamics
{

//--------------------------------------------------------------------------------------
// Base class for Transport Velocity Correction
//--------------------------------------------------------------------------------------
template <class BaseInteractionType>
class TransportVelocityCorrectionCKBase : public BaseInteractionType
{
  public:
    template <class DynamicsIdentifier>
    explicit TransportVelocityCorrectionCKBase(DynamicsIdentifier &identifier);
    virtual ~TransportVelocityCorrectionCKBase() {}

  protected:
    DiscreteVariable<Real> *dv_Vol_;                   ///< "VolumetricMeasure"
    DiscreteVariable<Vecd> *dv_dpos_;                  ///< "Position"
    DiscreteVariable<Vecd> *dv_zero_gradient_residue_; ///< "ZeroGradientResidue"
};

//--------------------------------------------------------------------------------------
// Main template declaration for TransportVelocityCorrectionCK
//--------------------------------------------------------------------------------------
template <typename...>
class TransportVelocityCorrectionCK;

template <class UpdatePolicy, class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
class TransportVelocityCorrectionCK<Inner<UpdatePolicy, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>
    : public TransportVelocityCorrectionCKBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = TransportVelocityCorrectionCKBase<Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;

  public:
    explicit TransportVelocityCorrectionCK(Inner<Parameters...> &inner_relation, Real coefficient = 0.2);

    virtual ~TransportVelocityCorrectionCK() {}

    //====================== Interact Kernel ======================//
    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Real *Vol_;
        Vecd *dpos_;
        Vecd *zero_gradient_residue_;
        ParticleScopeTypeKernel within_scope_;
    };

    //====================== Update Kernel ======================//
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;

        Real correction_scaling_;
        ResolutionType h_ratio_;
        LimiterType limiter_;
        Vecd *dpos_;
        Vecd *zero_gradient_residue_;
        ParticleScopeTypeKernel within_scope_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    Real h_ref_;              ///< e.g. reference smoothing length
    Real correction_scaling_; ///< typically coefficient * h_ref^2
    ResolutionType h_ratio_;  ///< e.g. for adaptive resolution
    LimiterType limiter_;     ///< e.g. a limiter on the final correction step
    ParticleScopeTypeCK<ParticleScopeType> within_scope_method_;
};
//----------------------------------------------
//  2) Partial specialization for Contact<...>
//----------------------------------------------
template <class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
class TransportVelocityCorrectionCK<Contact<Wall, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>
    : public TransportVelocityCorrectionCKBase<Interaction<Contact<Parameters...>>>, public Interaction<Wall>
{
    using BaseInteraction = TransportVelocityCorrectionCKBase<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;

  public:
    explicit TransportVelocityCorrectionCK(Contact<Parameters...> &contact_relation);
    virtual ~TransportVelocityCorrectionCK() {}

    //====================== Interact Kernel ======================//
    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Vecd *zero_gradient_residue_;
        Real *contact_wall_Vol_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    StdVec<DiscreteVariable<Real> *> dv_contact_wall_Vol_;
};

//--------------------------------------------------------------------------------------
// Alias Definitions for Specific Configurations
//--------------------------------------------------------------------------------------

using TransportVelocityCorrectionInnerNoCorrectionBulkParticlesCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, NoKernelCorrectionCK, SingleResolution, NoLimiter, BulkParticles>>;

using TransportVelocityCorrectionWallNoCorrectionBulkParticlesCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, NoKernelCorrectionCK, SingleResolution, NoLimiter, BulkParticles>,
        Contact<Wall, NoKernelCorrectionCK, SingleResolution, NoLimiter, BulkParticles>>;

using TransportVelocityLimitedCorrectionCorrectedComplexBulkParticlesCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, LinearCorrectionCK, SingleResolution, TruncatedLinear, BulkParticles>,
        Contact<Wall, LinearCorrectionCK, SingleResolution, TruncatedLinear, BulkParticles>>;

using TransportVelocityLimitedCorrectionCorrectedComplexBulkParticlesCKWithoutUpdate =
    TransportVelocityCorrectionCK<
        Inner<Base, LinearCorrectionCK, SingleResolution, TruncatedLinear, AllParticles>,
        Contact<Wall, LinearCorrectionCK, SingleResolution, TruncatedLinear, AllParticles>>;
using TransportVelocityCorrectedComplexBulkParticlesCKWithoutUpdate =
    TransportVelocityCorrectionCK<
        Inner<Base, NoKernelCorrectionCK, SingleResolution, NoLimiter, AllParticles>,
        Contact<Wall, NoKernelCorrectionCK, SingleResolution, NoLimiter, AllParticles>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // TRANSPORT_VELOCITY_CORRECTION_CK_H
