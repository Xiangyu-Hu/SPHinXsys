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
 * @file surface indication.h
 * @brief 	Here, we define the algorithm classes for indicating fluid surfaces.
 * @details Fluid indicators are mainly used here to classify different region in a fluid.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SURFACE_INDICATION_H
#define SURFACE_INDICATION_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class FreeSurfaceIndication
 * @brief  indicate the particles near the free surface of a fluid body.
 * Note that, SPHinXsys does not require this function for simulating general free surface flow problems.
 * However, some other applications may use this function, such as transport velocity formulation,
 * for masking some function which is only applicable for the bulk of the fluid body.
 */
template <class DataDelegationType>
class FreeSurfaceIndication : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit FreeSurfaceIndication(BaseRelationType &base_relation);
    virtual ~FreeSurfaceIndication(){};

  protected:
    StdLargeVec<int> &indicator_;
    StdLargeVec<Real> pos_div_;
    Real threshold_by_dimensions_;
};

/**
 * @class FreeSurfaceIndicationInner
 * @brief TBD.
 */
class FreeSurfaceIndicationInner : public FreeSurfaceIndication<FluidDataInner>
{
  public:
    explicit FreeSurfaceIndicationInner(BaseInnerRelation &inner_relation);
    virtual ~FreeSurfaceIndicationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real smoothing_length_;
    bool isVeryNearFreeSurface(size_t index_i);
};

/**
 * @class FreeSurfaceIndicationContact
 * @brief indicate the particles near the free fluid surface.
 */
class FreeSurfaceIndicationContact : public FreeSurfaceIndication<FluidContactData>
{
  public:
    FreeSurfaceIndicationContact(BaseContactRelation &contact_relation)
        : FreeSurfaceIndication<FluidContactData>(contact_relation){};
    virtual ~FreeSurfaceIndicationContact(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class SpatialTemporalFreeSurfaceIdentification
 * @brief using the spatial-temporal method to indicate the surface particles to avoid mis-judgement.
 */
template <class FreeSurfaceIdentificationType>
class SpatialTemporalFreeSurfaceIdentification : public FreeSurfaceIdentificationType
{
  public:
    template <typename... ConstructorArgs>
    explicit SpatialTemporalFreeSurfaceIdentification(ConstructorArgs &&...args);
    virtual ~SpatialTemporalFreeSurfaceIdentification(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<int> previous_surface_indicator_;
    bool isNearPreviousFreeSurface(size_t index_i);
};

/**
 * @class FreeSurfaceHeight
 * @brief Probe the free surface profile for a fluid body part by reduced operation.
 */
class FreeSurfaceHeight : public BaseLocalDynamicsReduce<Real, ReduceMax, BodyPartByCell>,
                          public FluidDataSimple
{
  protected:
    StdLargeVec<Vecd> &pos_;

  public:
    FreeSurfaceHeight(BodyPartByCell &body_part)
        : BaseLocalDynamicsReduce<Real, ReduceMax, BodyPartByCell>(body_part, Real(MinRealNumber)),
          FluidDataSimple(sph_body_), pos_(particles_->pos_)
    {
        quantity_name_ = "FreeSurfaceHeight";
    }
    virtual ~FreeSurfaceHeight(){};

    Real reduce(size_t index_i, Real dt = 0.0) { return pos_[index_i][1]; };
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // SURFACE_INDICATION_H
