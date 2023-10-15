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
 * @file general_constraint.h
 * @brief Particles are constrained on their position according to
 * different criteria.
 * @author	Xiangyu Hu
 */

#ifndef GENERAL_CONSTRAINT_H
#define GENERAL_CONSTRAINT_H

#include "base_general_dynamics.h"

namespace SPH
{
class LevelSetShape;

/**
 * @class ShapeSurfaceBounding
 * @brief constrain surface particles by
 * map constrained particles to geometry face and
 * r = r + phi * norm (vector distance to face)
 */
class ShapeSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>,
                             public GeneralDataDelegateSimple
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
} // namespace SPH
#endif // GENERAL_CONSTRAINT_H
