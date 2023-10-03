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
 * @file density_summation.h
 * @brief Here, we define the algorithm classes for computing
 * the density of a continuum by kernel function summation.
 * @details TBD.
 * @author Xiangyu Hu
 */

#ifndef DENSITY_SUMMATION_INNER_H
#define DENSITY_SUMMATION_INNER_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
class Inner;
class InnerAdaptive;
class Contact;
class ContactAdaptive;
template <typename InteractionType>
class FreeSurface;
template <typename InteractionType>
class FreeStream;

/**
 * @class ComplexInteraction
 * @brief A class that integrates multiple local dynamics.
 * Typically, it includes an inner interaction and one or
 * several contact interaction ad boundary conditions.
 */
template <typename... InteractionTypes>
class NewComplexInteraction;

template <template <typename... InteractionType> class LocalDynamicsName>
class NewComplexInteraction<LocalDynamicsName<>>
{
  public:
    NewComplexInteraction(){};

    void interaction(size_t index_i, Real dt = 0.0){};
};

template <template <typename... InteractionType> class LocalDynamicsName, class FirstInteraction, class... OtherInteractions>
class NewComplexInteraction<LocalDynamicsName<FirstInteraction, OtherInteractions...>> : public LocalDynamicsName<FirstInteraction>
{
  protected:
    NewComplexInteraction<LocalDynamicsName<OtherInteractions...>> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit NewComplexInteraction(FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets)
        : LocalDynamicsName<FirstInteraction>(first_parameter_set),
          other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        LocalDynamicsName<FirstInteraction>::interaction(index_i, dt);
        other_interactions_.interaction(index_i, dt);
    };
};

template <class DataDelegationType>
class BaseDensitySummation : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseDensitySummation(BaseRelationType &base_relation);
    virtual ~BaseDensitySummation(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_, &rho_sum_;
    Real rho0_, inv_sigma0_, W0_;
};

class DensitySummationInnerCommon : public BaseDensitySummation<FluidDataInner>
{
  public:
    explicit DensitySummationInnerCommon(BaseInnerRelation &inner_relation)
        : BaseDensitySummation<FluidDataInner>(inner_relation){};
    virtual ~DensitySummationInnerCommon(){};

  protected:
    void assignDensity(size_t index_i) { rho_[index_i] = rho_sum_[index_i]; };
    void reinitializeDensity(size_t index_i) { rho_[index_i] = SMAX(rho_sum_[index_i], rho0_); };
};

class DensitySummationContactCommon : public BaseDensitySummation<FluidContactData>
{
  public:
    explicit DensitySummationContactCommon(BaseContactRelation &contact_relation);
    virtual ~DensitySummationContactCommon(){};

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    Real ContactSummation(size_t index_i);
};

template <typename... InteractionTypes>
class DensitySummation;

template <>
class DensitySummation<Inner> : public DensitySummationInnerCommon
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation)
        : DensitySummationInnerCommon(inner_relation){};
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0) { assignDensity(index_i); };
};

template <>
class DensitySummation<InnerAdaptive> : public DensitySummationInnerCommon
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation);
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0) { assignDensity(index_i); };

  protected:
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    StdLargeVec<Real> &h_ratio_;
};

template <>
class DensitySummation<Contact> : public DensitySummationContactCommon
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation);
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <>
class DensitySummation<ContactAdaptive> : public DensitySummationContactCommon
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation);
    virtual ~DensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    StdLargeVec<Real> &h_ratio_;
};

template <class InteractionType>
class DensitySummation<FreeSurface<InteractionType>> : public DensitySummation<InteractionType>
{
  public:
    template <typename... ConstructorArgs>
    explicit DensitySummation(ConstructorArgs &&...args)
        : DensitySummation<InteractionType>(std::forward<ConstructorArgs>(args)...){};
    virtual ~DensitySummation(){};
    void update(size_t index_i, Real dt = 0.0)
    {
        DensitySummation<InteractionType>::reinitializeDensity(index_i);
    };
};

template <class InteractionType>
class DensitySummation<FreeStream<InteractionType>> : public DensitySummation<InteractionType>
{
  public:
    template <typename... ConstructorArgs>
    explicit DensitySummation(ConstructorArgs &&...args);
    virtual ~DensitySummation(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<int> &indicator_;
    bool isNearFreeSurface(size_t index_i);
};

template <class InnerInteractionType, class ContactInteractionType>
class BaseDensitySummationComplex
    : public NewComplexInteraction<DensitySummation<InnerInteractionType, ContactInteractionType>>
{
  public:
    explicit BaseDensitySummationComplex(ComplexRelation &complex_relation)
        : NewComplexInteraction<DensitySummation<InnerInteractionType, ContactInteractionType>>(
              complex_relation.getInnerRelation(), complex_relation.getContactRelation()){};
};
using DensitySummationComplex = BaseDensitySummationComplex<Inner, Contact>;
using DensitySummationComplexAdaptive = BaseDensitySummationComplex<InnerAdaptive, ContactAdaptive>;
using DensitySummationFreeSurfaceComplex = BaseDensitySummationComplex<FreeSurface<Inner>, Contact>;
using DensitySummationFreeSurfaceComplexAdaptive = BaseDensitySummationComplex<FreeSurface<InnerAdaptive>, ContactAdaptive>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_INNER_H