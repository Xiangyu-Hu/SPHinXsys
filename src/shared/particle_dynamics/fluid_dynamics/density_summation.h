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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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
 * @details We are using templates and their explicit or partial specializations
 * to identify variations of the interaction types..
 * @author Xiangyu Hu
 */

#ifndef DENSITY_SUMMATION_INNER_H
#define DENSITY_SUMMATION_INNER_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class DensitySummation;

template <class DataDelegationType>
class DensitySummation<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit DensitySummation(BaseRelationType &base_relation);
    virtual ~DensitySummation() {};

  protected:
    Real *rho_, *mass_, *rho_sum_, *Vol_;
    Real rho0_, inv_sigma0_, W0_;
};

template <>
class DensitySummation<Inner<Base>> : public DensitySummation<Base, DataDelegateInner>
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation)
        : DensitySummation<Base, DataDelegateInner>(inner_relation) {};
    virtual ~DensitySummation() {};
};

template <>
class DensitySummation<Inner<>> : public DensitySummation<Inner<Base>>
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation)
        : DensitySummation<Inner<Base>>(inner_relation) {};
    virtual ~DensitySummation() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using DensitySummationInner = DensitySummation<Inner<>>;

template <>
class DensitySummation<Inner<AdaptiveSmoothingLength>> : public DensitySummation<Inner<Base>>
{
  public:
    explicit DensitySummation(BaseInnerRelation &inner_relation);
    virtual ~DensitySummation() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    Real *h_ratio_;
};

template <>
class DensitySummation<Contact<Base>> : public DensitySummation<Base, DataDelegateContact>
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation);
    virtual ~DensitySummation() {};

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<Real *> contact_mass_;
    Real ContactSummation(size_t index_i);
};

template <>
class DensitySummation<Contact<>> : public DensitySummation<Contact<Base>>
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation)
        : DensitySummation<Contact<Base>>(contact_relation) {};
    virtual ~DensitySummation() {};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <>
class DensitySummation<Contact<AdaptiveSmoothingLength>> : public DensitySummation<Contact<Base>>
{
  public:
    explicit DensitySummation(BaseContactRelation &contact_relation);
    virtual ~DensitySummation() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    Real *h_ratio_;
};

template <typename... SummationType>
class DensitySummation<Inner<FreeSurface, SummationType...>> : public DensitySummation<Inner<SummationType...>>
{
  public:
    template <typename... Args>
    explicit DensitySummation(Args &&...args);
    virtual ~DensitySummation() {};
    void update(size_t index_i, Real dt = 0.0);
};
using DensitySummationFreeSurfaceInner = DensitySummation<Inner<FreeSurface>>;

struct FreeStream
{
    Real operator()(Real rho_sum, Real rho0, Real rho)
    {
        if (rho_sum < rho)
        {
            return rho_sum + SMAX(Real(0), (rho - rho_sum)) * rho0 / rho;
        }
        return rho_sum;
    };
};

struct NotNearSurface
{
    Real operator()(Real rho_sum, Real rho0, Real rho)
    {
        return rho;
    };
};

template <typename NearSurfaceType, typename... SummationType>
class DensitySummation<Inner<NearSurfaceType, SummationType...>>
    : public DensitySummation<Inner<SummationType...>>
{
  public:
    template <typename... Args>
    explicit DensitySummation(Args &&...args);
    virtual ~DensitySummation() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    NearSurfaceType near_surface_rho_;
    int *indicator_;
    bool isNearFreeSurface(size_t index_i);
};
using DensitySummationInnerNotNearSurface = DensitySummation<Inner<NotNearSurface>>;
using DensitySummationInnerFreeStream = DensitySummation<Inner<FreeStream>>;

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseDensitySummationComplex = ComplexInteraction<DensitySummation<InnerInteractionType, ContactInteractionTypes...>>;

using DensitySummationComplex = BaseDensitySummationComplex<Inner<>, Contact<>>;
using DensitySummationComplexAdaptive = BaseDensitySummationComplex<Inner<AdaptiveSmoothingLength>, Contact<AdaptiveSmoothingLength>>;
using DensitySummationComplexFreeSurface = BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>>;
using DensitySummationFreeSurfaceComplexAdaptive = BaseDensitySummationComplex<Inner<FreeSurface, AdaptiveSmoothingLength>, Contact<AdaptiveSmoothingLength>>;
using DensitySummationFreeStreamComplex = BaseDensitySummationComplex<Inner<FreeStream>, Contact<>>;
using DensitySummationFreeStreamComplexAdaptive = BaseDensitySummationComplex<Inner<FreeStream, AdaptiveSmoothingLength>, Contact<AdaptiveSmoothingLength>>;
using DensitySummationNotNearSurfaceComplex = BaseDensitySummationComplex<Inner<NotNearSurface>, Contact<>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_INNER_H