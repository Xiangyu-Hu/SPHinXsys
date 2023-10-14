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
 * @file 	relax_stepping.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Xiangyu Hu
 */

#ifndef RELAX_STEPPING_H
#define RELAX_STEPPING_H

#include "base_relax_dynamics.h"

namespace SPH
{
class GeometryShape;
class LevelSetShape;

namespace relax_dynamics
{
/**
 * @class GetTimeStepSizeSquare
 * @brief relaxation dynamics for particle initialization
 * computing the square of time step size
 */
class GetTimeStepSizeSquare : public LocalDynamicsReduce<Real, ReduceMax>,
                              public RelaxDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &acc_;
    Real h_ref_;

  public:
    explicit GetTimeStepSizeSquare(SPHBody &sph_body);
    virtual ~GetTimeStepSizeSquare(){};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value);
};

template <typename... InteractionTypes>
class ParticleRelaxation;

template <class DataDelegationType>
class ParticleRelaxation<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit ParticleRelaxation(BaseRelationType &base_relation);
    virtual ~ParticleRelaxation(){};

  protected:
    SPHAdaptation *sph_adaptation_;
    StdLargeVec<Vecd> &pos_;
};

template <>
class ParticleRelaxation<Inner<>>
    : public ParticleRelaxation<Base, RelaxDataDelegateInner>
{
  public:
    explicit ParticleRelaxation(BaseInnerRelation &inner_relation)
        : ParticleRelaxation<Base, RelaxDataDelegateInner>(inner_relation){};
    virtual ~ParticleRelaxation(){};
    void interaction(size_t index_i, Real dt_square);
};

template <>
class ParticleRelaxation<Inner<LevelSetCorrection>> : public ParticleRelaxation<Inner<>>
{
  public:
    explicit ParticleRelaxation(BaseInnerRelation &inner_relation);
    virtual ~ParticleRelaxation(){};
    void interaction(size_t index_i, Real dt_square);

  protected:
    LevelSetShape *level_set_shape_;
};

template <>
class ParticleRelaxation<Contact<>>
    : public ParticleRelaxation<Base, RelaxDataDelegateContact>
{
  public:
    explicit ParticleRelaxation(BaseContactRelation &contact_relation);
    virtual ~ParticleRelaxation(){};
    void interaction(size_t index_i, Real dt_square);
};

class UpdateSmoothingLengthRatioByShape : public LocalDynamics,
                                          public RelaxDataDelegateSimple
{
  protected:
    StdLargeVec<Real> &h_ratio_, &Vol_;
    StdLargeVec<Vecd> &pos_;
    Shape &target_shape_;
    ParticleRefinementByShape *particle_adaptation_;
    Real reference_spacing_;

  public:
    UpdateSmoothingLengthRatioByShape(SPHBody &sph_body, Shape &target_shape);
    explicit UpdateSmoothingLengthRatioByShape(SPHBody &sph_body);
    virtual ~UpdateSmoothingLengthRatioByShape(){};
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class ShapeSurfaceBounding
 * @brief constrain surface particles by
 * map constrained particles to geometry face and
 * r = r + phi * norm (vector distance to face)
 */
class ShapeSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>,
                             public RelaxDataDelegateSimple
{
  public:
    ShapeSurfaceBounding(NearShapeSurface &body_part);
    virtual ~ShapeSurfaceBounding(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    LevelSetShape *level_set_shape_;
    Real constrained_distance_;
};
} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_STEPPING_H
