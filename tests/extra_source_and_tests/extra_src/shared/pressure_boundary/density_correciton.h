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
 * @file 	density_correction.h
 * @brief 	Here, we define the algorithm of density correction for pressure boundary condition.
 * @author	Shuoguo Zhang and Xiangyu Hu
 */

#ifndef DENSITY_CORRECTION_H
#define DENSITY_CORRECTION_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class DensitySummationPressure;

template <class DataDelegationType>
class DensitySummationPressure<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit DensitySummationPressure(BaseRelationType &base_relation);
    virtual ~DensitySummationPressure(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_, &rho_sum_;
    Real rho0_, inv_sigma0_, W0_;
};

template <>
class DensitySummationPressure<Inner<Base>> : public DensitySummationPressure<Base, DataDelegateInner>
{
  public:
    explicit DensitySummationPressure(BaseInnerRelation &inner_relation)
        : DensitySummationPressure<Base, DataDelegateInner>(inner_relation){};
    virtual ~DensitySummationPressure(){};

  protected:
    void assignDensity(size_t index_i) { rho_[index_i] = rho_sum_[index_i]; };
};

template <>
class DensitySummationPressure<Inner<>> : public DensitySummationPressure<Inner<Base>>
{
  public:
    explicit DensitySummationPressure(BaseInnerRelation &inner_relation)
        : DensitySummationPressure<Inner<Base>>(inner_relation),
          buffer_particle_indicator_(*particles_->getVariableDataByName<int>("BufferParticleIndicator")){};
    virtual ~DensitySummationPressure(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0)
    {
        if (buffer_particle_indicator_[index_i] == 0)
            assignDensity(index_i);
    };

  protected:
    StdLargeVec<int> &buffer_particle_indicator_;
};

template <>
class DensitySummationPressure<Contact<Base>> : public DensitySummationPressure<Base, DataDelegateContact>
{
  public:
    explicit DensitySummationPressure(BaseContactRelation &contact_relation);
    virtual ~DensitySummationPressure(){};

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    Real ContactSummation(size_t index_i);
};

template <>
class DensitySummationPressure<Contact<>> : public DensitySummationPressure<Contact<Base>>
{
  public:
    explicit DensitySummationPressure(BaseContactRelation &contact_relation)
        : DensitySummationPressure<Contact<Base>>(contact_relation){};
    virtual ~DensitySummationPressure(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseDensitySummationPressureComplex = ComplexInteraction<DensitySummationPressure<InnerInteractionType, ContactInteractionTypes...>>;

using DensitySummationPressureComplex = BaseDensitySummationPressureComplex<Inner<>, Contact<>>;

} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_INNER_H