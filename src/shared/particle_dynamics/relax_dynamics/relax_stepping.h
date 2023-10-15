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
 * @file relax_stepping.h
 * @brief Particle relax is the process to produce body fitted
 * particle distribution with zero-order consistency.
 * @author	Xiangyu Hu
 */

#ifndef RELAX_STEPPING_H
#define RELAX_STEPPING_H

#include "base_relax_dynamics.h"
#include "general_constraint.h"

namespace SPH
{
class GeometryShape;
class LevelSetShape;

namespace relax_dynamics
{
template <typename... InteractionTypes>
class RelaxationResidue;

template <class DataDelegationType>
class RelaxationResidue<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit RelaxationResidue(BaseRelationType &base_relation);
    virtual ~RelaxationResidue(){};

  protected:
    SPHAdaptation *sph_adaptation_;
    StdLargeVec<Vecd> &residue_;
};

template <>
class RelaxationResidue<Inner<>>
    : public RelaxationResidue<Base, RelaxDataDelegateInner>
{
  public:
    explicit RelaxationResidue(BaseInnerRelation &inner_relation)
        : RelaxationResidue<Base, RelaxDataDelegateInner>(inner_relation){};
    virtual ~RelaxationResidue(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <>
class RelaxationResidue<Inner<LevelSetCorrection>> : public RelaxationResidue<Inner<>>
{
  public:
    explicit RelaxationResidue(BaseInnerRelation &inner_relation);
    virtual ~RelaxationResidue(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    LevelSetShape *level_set_shape_;
};

template <>
class RelaxationResidue<Contact<>>
    : public RelaxationResidue<Base, RelaxDataDelegateContact>
{
  public:
    explicit RelaxationResidue(BaseContactRelation &contact_relation);
    virtual ~RelaxationResidue(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class RelaxationScaling
 * @brief Obtain the scale for a particle relaxation step
 */
class RelaxationScaling : public LocalDynamicsReduce<Real, ReduceMax>,
                          public RelaxDataDelegateSimple
{
  public:
    explicit RelaxationScaling(SPHBody &sph_body);
    virtual ~RelaxationScaling(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value);

  protected:
    StdLargeVec<Vecd> &residue_;
    Real h_ref_;
};

/**
 * @class PositionRelaxation
 * @brief update the particle position for a relaxation step
 */
class PositionRelaxation : public LocalDynamics,
                           public RelaxDataDelegateSimple
{
  protected:
    SPHAdaptation *sph_adaptation_;
    StdLargeVec<Vecd> &pos_, &residue_;

  public:
    explicit PositionRelaxation(SPHBody &sph_body);
    virtual ~PositionRelaxation(){};
    void update(size_t index_i, Real scaling);
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

template <class RelaxationResidueType>
class RelaxationStep : public BaseDynamics<void>
{
  public:
    template <typename FirstArg, typename... OtherArgs>
    explicit RelaxationStep(FirstArg &first_arg, OtherArgs &&...other_args);
    virtual ~RelaxationStep(){};
    SimpleDynamics<ShapeSurfaceBounding> &SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

  protected:
    RealBody &real_body_;
    StdVec<SPHRelation *> &body_relations_;
    NearShapeSurface near_shape_surface_;
    InteractionDynamics<RelaxationResidueType> relaxation_residue_;
    ReduceDynamics<RelaxationScaling> relaxation_scaling_;
    SimpleDynamics<PositionRelaxation> position_relaxation_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

using RelaxationStepLevelSetCorrectionInner = RelaxationStep<RelaxationResidue<Inner<LevelSetCorrection>>>;
} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_STEPPING_H
