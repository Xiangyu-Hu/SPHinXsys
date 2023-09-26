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
 * @file 	surface_tension.h
 * @brief 	Here, we define the algorithm classes for fluid surface tensions.
 * @details Fluid indicators are mainly used here to classify different region in a fluid.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SURFACE_TENSION_H
#define SURFACE_TENSION_H

#include "base_fluid_dynamics.h"
#include "surface_indication.hpp"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class ColorFunctionGradient
 * @brief  indicate the particles near the interface of a fluid-fluid interaction and computing norm
 */

template <class DataDelegationType>
class ColorFunctionGradient : public FreeSurfaceIndication<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit ColorFunctionGradient(BaseRelationType &base_relation);
    virtual ~ColorFunctionGradient(){};

  protected:
    StdLargeVec<Vecd> color_grad_;
    StdLargeVec<Vecd> surface_norm_;
};

/**
 * @class ColorFunctionGradientInner
 * @brief  indicate the particles near the interface of a fluid-fluid interaction and computing norm
 */
class ColorFunctionGradientInner : public ColorFunctionGradient<FluidDataInner>
{
  public:
    explicit ColorFunctionGradientInner(BaseInnerRelation &inner_relation)
        : ColorFunctionGradient<FluidDataInner>(inner_relation){};
    virtual ~ColorFunctionGradientInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class ColorFunctionGradientInterpolationInner
 * @brief
 */
class ColorFunctionGradientInterpolationInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit ColorFunctionGradientInterpolationInner(BaseInnerRelation &inner_relation);
    virtual ~ColorFunctionGradientInterpolationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &Vol_;
    StdLargeVec<int> &indicator_;
    StdLargeVec<Vecd> &color_grad_;
    StdLargeVec<Vecd> &surface_norm_;
    StdLargeVec<Real> &pos_div_;
    Real threshold_by_dimensions_;
};

/**
 * @class SurfaceTensionAccelerationInner
 * @brief  the surface force induced acceleration
 */
class SurfaceTensionAccelerationInner : public LocalDynamics, public FluidDataInner
{
  public:
    SurfaceTensionAccelerationInner(BaseInnerRelation &inner_relation, Real gamma);
    explicit SurfaceTensionAccelerationInner(BaseInnerRelation &inner_relation);
    virtual ~SurfaceTensionAccelerationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real gamma_;
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<Vecd> &acc_prior_;
    StdLargeVec<int> &indicator_;
    StdLargeVec<Vecd> &color_grad_;
    StdLargeVec<Vecd> &surface_norm_;
};

/**
 * @class ColorFunctionGradientContact
 * @brief indicate the particles near the free fluid surface.
 */
class ColorFunctionGradientContact : public ColorFunctionGradient<FluidContactData>
{
  public:
    ColorFunctionGradientContact(BaseContactRelation &contact_relation)
        : ColorFunctionGradient<FluidContactData>(contact_relation){};
    virtual ~ColorFunctionGradientContact(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class 	SurfaceNormWithWall
 * @brief  Modify surface norm when contact with wall
 */
class SurfaceNormWithWall : public LocalDynamics, public FSIContactData
{
  public:
    SurfaceNormWithWall(BaseContactRelation &contact_relation, Real contact_angle);
    virtual ~SurfaceNormWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real contact_angle_;
    Real smoothing_length_;
    Real particle_spacing_;
    StdLargeVec<int> &indicator_;
    StdLargeVec<Vecd> &surface_norm_;
    StdLargeVec<Real> &pos_div_;
    StdVec<StdLargeVec<Vecd> *> wall_n_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // SURFACE_TENSION_H
