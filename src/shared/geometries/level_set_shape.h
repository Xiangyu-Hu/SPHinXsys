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
 * @file 	level_set_shape.h
 * @brief 	Here, we define geometry based on level set technique.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef LEVEL_SET_SHAPE_H
#define LEVEL_SET_SHAPE_H

#include "base_geometry.h"
#include "level_set.h"

#include <string>

namespace SPH
{
class IOEnvironment;
class SPHBody;
/**
 * @class LevelSetShape
 * @brief A shape using level set to define geometry
 */
class LevelSetShape : public Shape
{
  private:
    UniquePtrKeeper<BaseLevelSet> level_set_keeper_;
    SharedPtr<SPHAdaptation> sph_adaptation_;

  public:
    /** refinement_ratio is between body reference resolution and level set resolution */
    LevelSetShape(Shape &shape, SharedPtr<SPHAdaptation> sph_adaptation, Real refinement_ratio = 1.0);
    LevelSetShape(SPHBody &sph_body, Shape &shape, Real refinement_ratio = 1.0);

    virtual ~LevelSetShape(){};

    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override;
    virtual Vecd findClosestPoint(const Vecd &probe_point) override;

    Vecd findLevelSetGradient(const Vecd &probe_point);
    Real computeKernelIntegral(const Vecd &probe_point, Real h_ratio = 1.0);
    Vecd computeKernelGradientIntegral(const Vecd &probe_point, Real h_ratio = 1.0);
    /** small_shift_factor = 1.0 by default, can be increased for difficult geometries for smoothing */
    LevelSetShape *cleanLevelSet(Real small_shift_factor = 1.0);
    /** required to build level set from triangular mesh in stl file format. */
    LevelSetShape *correctLevelSetSign(Real small_shift_factor = 1.0);
    void writeLevelSet(IOEnvironment &io_environment);

  protected:
    BaseLevelSet &level_set_; /**< narrow bounded level set mesh. */

    virtual BoundingBox findBounds() override;
};
} // namespace SPH
#endif // LEVEL_SET_SHAPE_H
