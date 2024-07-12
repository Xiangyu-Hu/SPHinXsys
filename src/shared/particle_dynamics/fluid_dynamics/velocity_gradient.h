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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file velocity_gradient.h
 * @brief TBD.
 * @details TBD.
 * @author Xiangyu Hu
 */

#ifndef VELOCITY_GRADIENT_H
#define VELOCITY_GRADIENT_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class VelocityGradient;

template <class DataDelegationType>
class VelocityGradient<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit VelocityGradient(BaseRelationType &base_relation);
    virtual ~VelocityGradient(){};

  protected:
    Real *Vol_;
    Vecd *vel_;
    Matd *vel_grad_;
};

template <class KernelCorrectionType>
class VelocityGradient<Inner<KernelCorrectionType>>
    : public VelocityGradient<DataDelegateInner>
{
  public:
    explicit VelocityGradient(BaseInnerRelation &inner_relation);
    virtual ~VelocityGradient(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType kernel_correction_;
};
using VelocityGradientInner = VelocityGradient<Inner<NoKernelCorrection>>;

template <>
class VelocityGradient<Contact<Wall>> : public InteractionWithWall<VelocityGradient>
{
  public:
    explicit VelocityGradient(BaseContactRelation &wall_contact_relation);
    virtual ~VelocityGradient(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *distance_from_wall_;
};

template <class KernelCorrectionType>
using VelocityGradientWithWall = ComplexInteraction<VelocityGradient<Inner<KernelCorrectionType>, Contact<Wall>>>;

} // namespace fluid_dynamics
} // namespace SPH
#endif // VELOCITY_GRADIENT_H
