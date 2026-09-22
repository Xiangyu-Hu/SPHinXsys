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

struct NearFreeStream
{
    Real operator()(Real rho_sum, Real rho0, Real rho)
    {
      return (rho_sum < rho) ? (rho_sum + (rho - rho_sum) * rho0 / rho) : rho_sum;
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
using DensitySummationInnerFreeStream = DensitySummation<Inner<NearFreeStream>>;

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseDensitySummationComplex = ComplexInteraction<DensitySummation<InnerInteractionType, ContactInteractionTypes...>>;

using DensitySummationComplex = BaseDensitySummationComplex<Inner<>, Contact<>>;
using DensitySummationComplexAdaptive = BaseDensitySummationComplex<Inner<AdaptiveSmoothingLength>, Contact<AdaptiveSmoothingLength>>;
using DensitySummationComplexFreeSurface = BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>>;
using DensitySummationFreeSurfaceComplexAdaptive = BaseDensitySummationComplex<Inner<FreeSurface, AdaptiveSmoothingLength>, Contact<AdaptiveSmoothingLength>>;
using DensitySummationFreeStreamComplex = BaseDensitySummationComplex<Inner<NearFreeStream>, Contact<>>;
using DensitySummationFreeStreamComplexAdaptive = BaseDensitySummationComplex<Inner<NearFreeStream, AdaptiveSmoothingLength>, Contact<AdaptiveSmoothingLength>>;
using DensitySummationNotNearSurfaceComplex = BaseDensitySummationComplex<Inner<NotNearSurface>, Contact<>>;

/**
 * @class ShepardDensityRegularizationWithWall
 * @brief Zeroth-order consistent density regularization for a fluid with solid walls.
 * @details The density is reconstructed from the current fluid density field with
 * a normalized SPH interpolant. Solid-wall particles complete the denominator
 * with their reference volume and use a zero-normal-gradient extension of the
 * target fluid density. Applying the same expression to every real fluid
 * particle avoids a discontinuous switch at the edge of the wall contact list.
 *
 * This is an opt-in, fixed-smoothing-length operator; it does not replace the
 * existing density summation defaults. Use it with InteractionWithUpdate so
 * that all interactions read the old density and volume before any update.
 * Example: InteractionWithUpdate<ShepardDensityRegularizationWithWall>
 *              regularize_density(fluid_inner_relation, fluid_wall_relation);
 * The reconstructed result is staged in ShepardDensity. DensitySummation is
 * neither computed nor modified: it denotes the raw geometric density sum used
 * by existing boundary consumers such as FreeStreamVelocityCorrection. Callers
 * using those consumers must compute that raw sum separately before applying
 * this regularization; the regularization is not a replacement for that field.
 *
 * The zero-normal-gradient wall extension is not a hydrostatic or accelerating
 * wall pressure condition. It does not supply an adaptive-resolution or a
 * free-surface density constraint. Selecting this operator changes the density
 * evolution, so an existing trajectory regression database is not validation
 * of the new method. The caller must validate the chosen physical application.
 */
class ShepardDensityRegularizationWithWall
    : public LocalDynamics, public DataDelegateInner, public DataDelegateContact
{
  public:
    ShepardDensityRegularizationWithWall(BaseInnerRelation &inner_relation,
                                         BaseContactRelation &wall_contact_relation);
    virtual ~ShepardDensityRegularizationWithWall() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *rho_, *mass_, *rho_regularized_, *Vol_;
    Real W0_;
    StdVec<Real> contact_inv_rho0_;
    StdVec<Real *> contact_mass_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_INNER_H
