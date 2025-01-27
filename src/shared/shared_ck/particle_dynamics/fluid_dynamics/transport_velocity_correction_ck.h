#ifndef TRANSPORT_VELOCITY_CORRECTION_CK_H
#define TRANSPORT_VELOCITY_CORRECTION_CK_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"

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
    DiscreteVariable<Real> *dv_Vol_;                     ///< "VolumetricMeasure"
    DiscreteVariable<Vecd> *dv_dpos_;                    ///< "Position"
    DiscreteVariable<Vecd> *dv_zero_gradient_residue_;  ///< "ZeroGradientResidue"

};

//--------------------------------------------------------------------------------------
// Main template declaration for TransportVelocityCorrectionCK
//--------------------------------------------------------------------------------------
template <typename...>
class TransportVelocityCorrectionCK;

template <class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
class TransportVelocityCorrectionCK<Inner<WithUpdate, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>
    : public TransportVelocityCorrectionCKBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = TransportVelocityCorrectionCKBase<Interaction<Inner<Parameters...>>>;

    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

public:
    /// Constructor: accepts `Relation<Inner<>>` and optional coefficient.
    explicit TransportVelocityCorrectionCK(Relation<Inner<Parameters...>> &inner_relation);

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
        Vecd *dpos_, *zero_gradient_residue_;
        ParticleScopeType within_scope_;

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
        Vecd *dpos_, *zero_gradient_residue_;
        ParticleScopeType within_scope_;
    };

protected:
    KernelCorrectionType kernel_correction_; ///< The user-chosen kernel correction policy

    Real h_ref_;               ///< e.g. reference smoothing length
    Real correction_scaling_;  ///< typically coefficient * h_ref^2

    ResolutionType h_ratio_;   ///< e.g. for adaptive resolution
    LimiterType limiter_;      ///< e.g. a limiter on the final correction step
    ParticleScopeType within_scope_;

};
//----------------------------------------------
//  2) Partial specialization for Contact<...> 
//----------------------------------------------
template <class KernelCorrectionType, class ResolutionType, class LimiterType, class ParticleScopeType, typename... Parameters>
class TransportVelocityCorrectionCK<Contact<Wall, KernelCorrectionType, ResolutionType, LimiterType, ParticleScopeType, Parameters...>>
    : public TransportVelocityCorrectionCKBase<Interaction<Contact<Wall, Parameters...>>>
{
        using BaseInteraction = TransportVelocityCorrectionCKBase<Interaction<Contact<Wall, Parameters...>>>;
        using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

    public:
        explicit TransportVelocityCorrectionCK(Relation<Contact<Parameters...>> &contact_relation);
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
           Vecd  *zero_gradient_residue_;
           Real* contact_wall_Vol_;
        };
          protected:

    KernelCorrectionType kernel_correction_;
        StdVec<DiscreteVariable<Real> *> dv_contact_wall_Vol_;

};

//--------------------------------------------------------------------------------------
// Alias Definitions for Specific Configurations
//--------------------------------------------------------------------------------------

using TransportVelocityCorrectionInnerNoCorrectionCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, NoKernelCorrectionCK, SingleResolution, NoLimiter, SPH::BulkParticles>>;

using TransportVelocityCorrectionWallNoCorrectionCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, NoKernelCorrectionCK, SingleResolution, NoLimiter, SPH::BulkParticles>,
        Contact<Wall, NoKernelCorrectionCK, SingleResolution, NoLimiter, SPH::BulkParticles>>;

using TransportVelocityLimitedCorrectionCorrectedComplexCK =
    TransportVelocityCorrectionCK<
        Inner<WithUpdate, LinearCorrectionCK, SingleResolution, TruncatedLinear, SPH::BulkParticles>,
        Contact<Wall, LinearCorrectionCK, SingleResolution, TruncatedLinear, SPH::BulkParticles>>;
// using TransportVelocityCorrectionInnerLinearCorrectionCK =
//     TransportVelocityCorrectionCK<
//         Inner<LinearCorrectionCK, SingleResolution, NoLimiter, ParticleScope>>;

}  // namespace fluid_dynamics
}  // namespace SPH

#endif  // TRANSPORT_VELOCITY_CORRECTION_CK_H
