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
 * @file 	level_set_shape.h
 * @brief 	Here, we define geometry based on level set technique.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef LEVEL_SET_SHAPE_H
#define LEVEL_SET_SHAPE_H

#include "base_geometry.h"
#include "level_set.hpp"

namespace SPH
{
class SPHBody;
class SPHSystem;
/**
 * @class LevelSetShape
 * @brief A shape using level set to define geometry
 */
class LevelSetShape : public Shape
{
  private:
    UniquePtrKeeper<LevelSet> level_set_keeper_;
    SharedPtr<SPHAdaptation> sph_adaptation_;

  public:
    /** refinement is between body reference resolution and level set resolution */
    LevelSetShape(SPHBody &sph_body, Shape &shape, Real refinement = 1.0,
                  UsageType usage_type = UsageType::Volumetric);

    template <class ExecutionPolicy>
    LevelSetShape(const ExecutionPolicy &ex_policy, SPHSystem &sph_system, const SPHAdaptation &sph_adaptation,
                  Shape &shape, Real refinement = 1.0, UsageType usage_type = UsageType::Volumetric)
        : LevelSetShape(sph_system, sph_adaptation, shape, refinement)
    {
        finishInitialization(ex_policy, usage_type);
    };

    virtual ~LevelSetShape() {};

    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override;
    virtual Vecd findClosestPoint(const Vecd &probe_point) override;
    virtual BoundingBoxd findBounds() override;

    template <class ExecutionPolicy>
    void finishInitialization(const ExecutionPolicy &ex_policy, UsageType usage_type)
    {
        level_set_.finishInitialization(ex_policy, usage_type);
    };
    Vecd findLevelSetGradient(const Vecd &probe_point);
    Real computeKernelIntegral(const Vecd &probe_point, Real h_ratio = 1.0);
    Vecd computeKernelGradientIntegral(const Vecd &probe_point, Real h_ratio = 1.0);
    Matd computeKernelSecondGradientIntegral(const Vecd &probe_point, Real h_ratio = 1.0);
    /** small_shift_factor = 1.0 by default, can be increased for difficult geometries for smoothing */
    LevelSetShape *cleanLevelSet(UnsignedInt repeat_times = 1);
    /** required to build level set from triangular mesh in stl file format. */
    LevelSetShape *correctLevelSetSign();
    LevelSetShape *writeLevelSet();
    LevelSet &getLevelSet() { return level_set_; }

    template <typename DataType>
    LevelSetShape *addPackageVariableToWrite(const std::string &variable_name)
    {
        level_set_.addPackageVariableToWrite<DataType>(variable_name);
        return this;
    };

    template <typename DataType>
    LevelSetShape *addCellVariableToWrite(const std::string &variable_name)
    {
        level_set_.addCellVariableToWrite<DataType>(variable_name);
        return this;
    };

  protected:
    LevelSetShape(SPHSystem &sph_system, const SPHAdaptation &sph_adaptation, Shape &shape, Real refinement);
    SPHSystem &sph_system_; /**< for write level to file. */
    LevelSet &level_set_;   /**< narrow bounded level set mesh. */
};
} // namespace SPH
#endif // LEVEL_SET_SHAPE_H
