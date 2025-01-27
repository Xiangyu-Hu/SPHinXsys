#ifndef TRANSPORT_VELOCITY_CORRECTION_CK_H
#define TRANSPORT_VELOCITY_CORRECTION_CK_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{
/*
 * ParticleScopeTypeCK determines the scope of particles for which transport 
 * velocity correction is applied. It is a template class specialized for 
 * different particle types, defining rules for inclusion in the computation.
 */
template <typename ScopeType>
class ParticleScopeTypeCK; 
//-------------------------------------------------------------------------------------------------
// 1) Specialization for AllParticles
//-------------------------------------------------------------------------------------------------
template <>
class ParticleScopeTypeCK<AllParticles> : public WithinScope
{
public:
    // Constructor
    explicit ParticleScopeTypeCK(BaseParticles* particles)
        : WithinScope() 
    {}

    // Nested functor-like class:
    class ComputingKernel
    {
    public:
        // The signature typically follows the style of other SPH "ComputingKernel" constructors
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleScopeTypeCK<AllParticles> &encloser,
                        ComputingKernelType & computing_kernel)
        {};

        bool operator()(size_t index_i) const
        {
            return true; // Always in scope
        };
    };
};

//-------------------------------------------------------------------------------------------------
// 2) Specialization for BulkParticles (which is typically IndicatedParticles<0> in SPH code)
//    i.e. we only want particles that have indicator_[i] == 0
//-------------------------------------------------------------------------------------------------
template <>
class ParticleScopeTypeCK<BulkParticles> : public WithinScope
{
public:
    explicit ParticleScopeTypeCK(BaseParticles* particles)
        : WithinScope(),
          dv_indicator_(particles->getVariableByName<int>("Indicator"))
    {}

    // Nested class implementing the boolean check
    class ComputingKernel
    {
    public:
        // Typically, we pass the "encloser" object to get the pointer from dv_indicator_
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleScopeTypeCK<BulkParticles> &encloser,
                        ComputingKernelType & computing_kernel)
            : indicator_(encloser.dv_indicator_->DelegatedData(ex_policy))
        {}

        bool operator()(size_t index_i) const
        {
            return (indicator_[index_i] == 0);
        }

    protected:
        int* indicator_;
    };

protected:
    DiscreteVariable<int> *dv_indicator_;
};

//-------------------------------------------------------------------------------------------------
// 3) Specialization for NotIndicatedParticles (which is typically "indicator != 0")
//-------------------------------------------------------------------------------------------------
template <>
class ParticleScopeTypeCK<NotIndicatedParticles<0>> : public WithinScope
{
public:
    explicit ParticleScopeTypeCK(BaseParticles* particles)
        : WithinScope(),
          dv_indicator_(particles->getVariableByName<int>("Indicator"))
    {
    }

    class ComputingKernel
    {
    public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        ParticleScopeTypeCK<NotIndicatedParticles<0>> &encloser,
                        ComputingKernelType & computing_kernel)
            : indicator_(encloser.dv_indicator_->DelegatedData(ex_policy))
        {}

        bool operator()(size_t index_i) const
        {
            return (indicator_[index_i] != 0);
        }

    protected:
        int* indicator_;
    };

protected:
    DiscreteVariable<int> *dv_indicator_;
};


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
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;



public:    
    explicit TransportVelocityCorrectionCK(Relation<Inner<Parameters...>> &inner_relation, Real coefficient = 0.2);
    
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
        Vecd *dpos_, *zero_gradient_residue_;
        ParticleScopeTypeKernel within_scope_;
    };

protected:
    KernelCorrectionType kernel_correction_;
    Real h_ref_;               ///< e.g. reference smoothing length
    Real correction_scaling_;  ///< typically coefficient * h_ref^2
    ResolutionType h_ratio_;   ///< e.g. for adaptive resolution
    LimiterType limiter_;      ///< e.g. a limiter on the final correction step
    ParticleScopeTypeCK<ParticleScopeType> within_scope_method_;

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
        using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;

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


}  // namespace fluid_dynamics
}  // namespace SPH

#endif  // TRANSPORT_VELOCITY_CORRECTION_CK_H
