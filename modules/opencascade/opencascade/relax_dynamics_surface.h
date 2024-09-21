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
 * @file 	relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef RELAX_DYNAMICS_SURFACE_H
#define RELAX_DYNAMICS_SURFACE_H

#include "sphinxsys.h"
#include "surface_shape.h"
#include "vector.h"

namespace SPH
{

class SurfaceShape;

namespace relax_dynamics
{
class ShapeSurfaceConstraint : public LocalDynamics
{
  public:
    ShapeSurfaceConstraint(RealBody &real_body, Shape &shape);
    virtual ~ShapeSurfaceConstraint(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *pos_;
    Shape *shape_;
};

class RelaxationStepInnerFirstHalf : public BaseDynamics<void>
{
  public:
    explicit RelaxationStepInnerFirstHalf(BaseInnerRelation &inner_relation, Shape &shape);
    virtual ~RelaxationStepInnerFirstHalf(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    RealBody *real_body_;
    BaseInnerRelation &inner_relation_;
    InteractionDynamics<RelaxationResidue<Inner<>>> relaxation_acceleration_inner_;
};

class RelaxationStepInnerSecondHalf : public BaseDynamics<void>
{
  public:
    explicit RelaxationStepInnerSecondHalf(BaseInnerRelation &inner_relation, Shape &shape);
    virtual ~RelaxationStepInnerSecondHalf(){};
    SimpleDynamics<ShapeSurfaceConstraint> &SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

  protected:
    RealBody *real_body_;
    ReduceDynamics<RelaxationScaling> get_time_step_square_;
    SimpleDynamics<PositionRelaxation> update_particle_position_;
    SimpleDynamics<ShapeSurfaceConstraint> surface_bounding_;
};

/**
 * @class SurfaceNormalDirection
 * @brief get the normal direction of surface particles.
 */
class SurfaceNormalDirection : public LocalDynamics
{
  public:
    explicit SurfaceNormalDirection(SPHBody &sph_body, Shape &shape);
    virtual ~SurfaceNormalDirection(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    SurfaceShape *surface_shape_;
    Vecd *pos_, *n_;
};

} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_DYNAMICS_H
