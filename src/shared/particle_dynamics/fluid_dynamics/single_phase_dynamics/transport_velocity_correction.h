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
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef TRANSPORT_VELOCITY_CORRECTION_H
#define TRANSPORT_VELOCITY_CORRECTION_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class TransportVelocityCorrection
 * @brief TBD
 * @details If single (acoustic) time step is used, the coefficient should be decrease
 * to about 1/4 of the default value.
 */
template <class DataDelegationType, class KernelCorrectionType, class ResolutionType, class ParticleScope>
class TransportVelocityCorrection : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit TransportVelocityCorrection(BaseRelationType &base_relation, Real coefficient);
    virtual ~TransportVelocityCorrection(){};

  protected:
    const Real correction_scaling_;
    StdLargeVec<Vecd> &pos_;
    KernelCorrectionType kernel_correction_;
    ResolutionType h_ratio_;
    ParticleScope checkWithinScope;
};

template <class KernelCorrectionType, class ResolutionType, class ParticleScope>
class TransportVelocityCorrectionInner
    : public TransportVelocityCorrection<FluidDataInner, KernelCorrectionType, ResolutionType, ParticleScope>
{
  public:
    explicit TransportVelocityCorrectionInner(BaseInnerRelation &inner_relation, Real coefficient = 0.2);
    TransportVelocityCorrectionInner(LocalDynamicsParameters<BaseInnerRelation, Real> parameters)
        : TransportVelocityCorrectionInner(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~TransportVelocityCorrectionInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <class KernelCorrectionType, class ResolutionType, class ParticleScope>
class TransportVelocityCorrectionWithBoundary
    : public TransportVelocityCorrection<FluidContactData, KernelCorrectionType, ResolutionType, ParticleScope>
{
  public:
    explicit TransportVelocityCorrectionWithBoundary(BaseContactRelation &contact_relation, Real coefficient = 0.2);
    TransportVelocityCorrectionWithBoundary(LocalDynamicsParameters<BaseContactRelation, Real> parameters)
        : TransportVelocityCorrectionWithBoundary(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~TransportVelocityCorrectionWithBoundary(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <class KernelCorrectionType, class ResolutionType, class ParticleScope>
class TransportVelocityCorrectionContact
    : public TransportVelocityCorrection<FluidContactData, KernelCorrectionType, ResolutionType, ParticleScope>
{
  public:
    explicit TransportVelocityCorrectionContact(BaseContactRelation &contact_relation, Real coefficient = 0.2);
    TransportVelocityCorrectionContact(LocalDynamicsParameters<BaseContactRelation, Real> parameters)
        : TransportVelocityCorrectionContact(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~TransportVelocityCorrectionContact(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<KernelCorrectionType> contact_kernel_corrections_;
};

template <class ParticleScope>
class TransportVelocityCorrectionComplex
    : public OldComplexInteraction<TransportVelocityCorrectionInner<NoKernelCorrection, SingleResolution, ParticleScope>,
                                TransportVelocityCorrectionWithBoundary<NoKernelCorrection, SingleResolution, ParticleScope>>
{
  public:
    explicit TransportVelocityCorrectionComplex(ComplexRelation &complex_relation, Real coefficient = 0.2)
        : OldComplexInteraction<TransportVelocityCorrectionInner<NoKernelCorrection, SingleResolution, ParticleScope>,
                             TransportVelocityCorrectionWithBoundary<NoKernelCorrection, SingleResolution, ParticleScope>>(
              LocalDynamicsParameters(complex_relation.getInnerRelation(), coefficient),
              LocalDynamicsParameters(complex_relation.getContactRelation(), coefficient)){};
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // TRANSPORT_VELOCITY_CORRECTION_H
