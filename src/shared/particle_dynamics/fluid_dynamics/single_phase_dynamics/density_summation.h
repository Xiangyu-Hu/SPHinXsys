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
/**
 * @class BaseDensitySummation
 * @brief Base class for computing density by summation
 */
template <class DataDelegationType>
class BaseDensitySummation : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseDensitySummation(BaseRelationType &base_relation);
    virtual ~BaseDensitySummation(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_, &rho_sum_;
    Real rho0_, inv_sigma0_;
};

/**
 * @class BaseDensitySummationInner
 * @brief Base class for computing density by summation
 */
class BaseDensitySummationInner : public BaseDensitySummation<FluidDataInner>
{
  public:
    explicit BaseDensitySummationInner(BaseInnerRelation &inner_relation)
        : BaseDensitySummation<FluidDataInner>(inner_relation){};
    virtual ~BaseDensitySummationInner(){};
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class DensitySummationInner
 * @brief  computing density by summation
 */
class DensitySummationInner : public BaseDensitySummationInner
{
  public:
    explicit DensitySummationInner(BaseInnerRelation &inner_relation);
    virtual ~DensitySummationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real W0_;
};

/**
 * @class DensitySummationInnerAdaptive
 * @brief  computing density by summation with variable smoothing length
 */
class DensitySummationInnerAdaptive : public BaseDensitySummationInner
{
  public:
    explicit DensitySummationInnerAdaptive(BaseInnerRelation &inner_relation);
    virtual ~DensitySummationInnerAdaptive(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class BaseDensitySummationContact
 * @brief computing density by summation considering contribution from contact bodies
 */
class BaseDensitySummationContact : public BaseDensitySummation<FluidContactData>
{
  public:
    explicit BaseDensitySummationContact(BaseContactRelation &contact_relation);
    virtual ~BaseDensitySummationContact(){};

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    Real ContactSummation(size_t index_i);
};

/**
 * @class DensitySummationContact
 * @brief computing density by summation considering contribution from contact bodies
 */
class DensitySummationContact : public BaseDensitySummationContact
{
  public:
    explicit DensitySummationContact(BaseContactRelation &contact_relation)
        : BaseDensitySummationContact(contact_relation){};
    virtual ~DensitySummationContact(){};

    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class DensitySummationContactAdaptive
 * @brief computing density by summation considering  contribution from contact bodies
 */
class DensitySummationContactAdaptive : public BaseDensitySummationContact
{
  public:
    explicit DensitySummationContactAdaptive(BaseContactRelation &contact_relation);
    virtual ~DensitySummationContactAdaptive(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class DensitySummationFreeSurface
 * @brief computing density by summation with a re-normalization for free surface flows
 */
template <class DensitySummationType>
class DensitySummationFreeSurface : public DensitySummationType
{
  public:
    template <typename... ConstructorArgs>
    explicit DensitySummationFreeSurface(ConstructorArgs &&...args)
        : DensitySummationType(std::forward<ConstructorArgs>(args)...){};
    virtual ~DensitySummationFreeSurface(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real ReinitializedDensity(Real rho_sum, Real rho_0)
    {
        return SMAX(rho_sum, rho_0);
    };
};

/**
 * @class DensitySummationFreeStream
 * @brief The density is smoothed if the particle is near fluid surface.
 */
template <class DensitySummationFreeSurfaceType>
class DensitySummationFreeStream : public DensitySummationFreeSurfaceType
{
  public:
    template <typename... ConstructorArgs>
    explicit DensitySummationFreeStream(ConstructorArgs &&...args);
    virtual ~DensitySummationFreeStream(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<int> &indicator_;
    bool isNearFreeSurface(size_t index_i);
};

class DensitySummationFreeSurfaceComplex
    : public DensitySummationFreeSurface<ComplexInteraction<DensitySummationInner, DensitySummationContact>>
{
  public:
    explicit DensitySummationFreeSurfaceComplex(ComplexRelation &fluid_wall_relation)
        : DensitySummationFreeSurface<ComplexInteraction<DensitySummationInner, DensitySummationContact>>(
              fluid_wall_relation.getInnerRelation(), fluid_wall_relation.getContactRelation()){};
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_INNER_H