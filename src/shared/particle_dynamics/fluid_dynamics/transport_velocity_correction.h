/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file transport_velocity_correction.h
 * @brief The particle positions are corrected for more uniformed distribution
 * when there is negative pressure in the flow.
 * @details Note that the default coefficient is for using the dual time criteria method:
 * Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 * C Zhang, M Rezavand, X Hu - Journal of Computational Physics,
 * Volume 404, 1 March 2020, 109135.
 * If single (acoustic) time step is used, the coefficient should be decrease
 * to about 1/4 of the default value.
 * @author	Xiangyu Hu
 */

#ifndef TRANSPORT_VELOCITY_CORRECTION_H
#define TRANSPORT_VELOCITY_CORRECTION_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... T>
class TransportVelocityCorrection;

template <class DataDelegationType, class KernelCorrectionType, class ParticleScope>
class TransportVelocityCorrection<Base, DataDelegationType, KernelCorrectionType, ParticleScope>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit TransportVelocityCorrection(BaseRelationType &base_relation);
    virtual ~TransportVelocityCorrection() {};

  protected:
    Vecd *zero_gradient_residue_;
    KernelCorrectionType kernel_correction_;
    ParticleScope within_scope_;
};

template <class ResolutionType, class LimiterType, typename... CommonControlTypes>
class TransportVelocityCorrection<Inner<ResolutionType, LimiterType>, CommonControlTypes...>
    : public TransportVelocityCorrection<Base, DataDelegateInner, CommonControlTypes...>
{
  public:
    explicit TransportVelocityCorrection(BaseInnerRelation &inner_relation, Real coefficient = 0.2);
    template <typename BodyRelationType, typename FirstArg>
    explicit TransportVelocityCorrection(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : TransportVelocityCorrection(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~TransportVelocityCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    const Real h_ref_, correction_scaling_;
    Real *Vol_;
    Vecd *pos_;
    ResolutionType h_ratio_;
    LimiterType limiter_;
};
template <class LimiterType, class ParticleScope>
using TransportVelocityCorrectionInner =
    TransportVelocityCorrection<Inner<SingleResolution, LimiterType>, NoKernelCorrection, ParticleScope>;

template <typename... CommonControlTypes>
class TransportVelocityCorrection<Contact<Boundary>, CommonControlTypes...>
    : public TransportVelocityCorrection<Base, DataDelegateContact, CommonControlTypes...>
{
  public:
    explicit TransportVelocityCorrection(BaseContactRelation &contact_relation);
    virtual ~TransportVelocityCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real *> wall_Vol_;
};

template <class KernelCorrectionType, typename... CommonControlTypes>
class TransportVelocityCorrection<Contact<>, KernelCorrectionType, CommonControlTypes...>
    : public TransportVelocityCorrection<Base, DataDelegateContact, KernelCorrectionType, CommonControlTypes...>
{
  public:
    explicit TransportVelocityCorrection(BaseContactRelation &contact_relation);
    virtual ~TransportVelocityCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<KernelCorrectionType> contact_kernel_corrections_;
    StdVec<Real *> contact_Vol_;
};

template <class ResolutionType, class LimiterType, typename... CommonControlTypes>
using BaseTransportVelocityCorrectionComplex =
    ComplexInteraction<TransportVelocityCorrection<Inner<ResolutionType, LimiterType>, Contact<Boundary>>, CommonControlTypes...>;

template <class ParticleScope>
using TransportVelocityCorrectionComplex =
    BaseTransportVelocityCorrectionComplex<SingleResolution, NoLimiter, NoKernelCorrection, ParticleScope>;

template <class ParticleScope>
using TransportVelocityCorrectionCorrectedComplex =
    BaseTransportVelocityCorrectionComplex<SingleResolution, NoLimiter, LinearGradientCorrection, ParticleScope>;

template <class ParticleScope>
using TransportVelocityCorrectionCorrectedForOpenBoundaryFlowComplex =
    BaseTransportVelocityCorrectionComplex<SingleResolution, NoLimiter, LinearGradientCorrectionWithBulkScope, ParticleScope>;

template <class ParticleScope>
using TransportVelocityLimitedCorrectionCorrectedForOpenBoundaryFlowComplex =
    BaseTransportVelocityCorrectionComplex<SingleResolution, TruncatedLinear, LinearGradientCorrectionWithBulkScope, ParticleScope>;

template <class ParticleScope>
using TransportVelocityLimitedCorrectionComplex =
    BaseTransportVelocityCorrectionComplex<SingleResolution, TruncatedLinear, NoKernelCorrection, ParticleScope>;

template <class ParticleScope>
using TransportVelocityLimitedCorrectionCorrectedComplex =
    BaseTransportVelocityCorrectionComplex<SingleResolution, TruncatedLinear, LinearGradientCorrection, ParticleScope>;

template <class ParticleScope>
using TransportVelocityCorrectionComplexAdaptive =
    BaseTransportVelocityCorrectionComplex<AdaptiveResolution, NoLimiter, NoKernelCorrection, ParticleScope>;

template <class ResolutionType, typename... CommonControlTypes>
using BaseMultiPhaseTransportVelocityCorrectionComplex =
    ComplexInteraction<TransportVelocityCorrection<Inner<ResolutionType, NoLimiter>, Contact<>, Contact<Boundary>>, CommonControlTypes...>;

template <class ParticleScope>
using MultiPhaseTransportVelocityCorrectionComplex =
    BaseMultiPhaseTransportVelocityCorrectionComplex<SingleResolution, NoKernelCorrection, ParticleScope>;

} // namespace fluid_dynamics
} // namespace SPH
#endif // TRANSPORT_VELOCITY_CORRECTION_H
