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
 * @file 	complex_shape.h
 * @brief 	Here, we define the a container of different type of shapes, which allow
 * 			the boolean operation between the different type of shapes. The shapes
 * 			can be defined previously and add to this complex shapes container.
 */

#ifndef COMPLEX_SHAPE_H
#define COMPLEX_SHAPE_H

#include "base_geometry.h"
#include "geometric_shape.h"
#include "level_set_shape.h"
#include "transform_shape.h"

namespace SPH
{
class LevelSetShape;

class ComplexShape : public BinaryShapes
{
  public:
    explicit ComplexShape(const std::string &shape_name)
        : BinaryShapes(shape_name){};
    virtual ~ComplexShape(){};

    template <typename... Args>
    LevelSetShape *defineLevelSetShape(SPHBody &sph_body, const std::string &shape_name, Args &&...args)
    {
        size_t index = getSubShapeIndexByName(shape_name);
        LevelSetShape *level_set_shape = sub_shape_ptrs_keeper_[index].createPtr<LevelSetShape>(
            sph_body, *sub_shapes_and_ops_[index].first, std::forward<Args>(args)...);
        sub_shapes_and_ops_[index].first = DynamicCast<Shape>(this, level_set_shape);
        return level_set_shape;
    };
};

using DefaultShape = ComplexShape;

/**
 * @class AlignedBoxShape
 * @brief Used to describe a bounding box in which
 * the upper bound direction is aligned to the normal of a planar piece on a shape.
 */
class AlignedBoxShape : public TransformShape<GeometricShapeBox>
{
    const int alignment_axis_;

  public:
    /** construct directly */
    template <typename... Args>
    explicit AlignedBoxShape(int upper_bound_axis, const Transform &transform, Args &&...args)
        : TransformShape<GeometricShapeBox>(transform, std::forward<Args>(args)...),
          alignment_axis_(upper_bound_axis){};
    /** construct from a shape already has aligned boundaries */
    template <typename... Args>
    explicit AlignedBoxShape(int upper_bound_axis, const Shape &shape, Args &&...args)
        : TransformShape<GeometricShapeBox>(
              Transform(Vecd(0.5 * (shape.bounding_box_.second_ + shape.bounding_box_.first_))),
              0.5 * (shape.bounding_box_.second_ - shape.bounding_box_.first_), std::forward<Args>(args)...),
          alignment_axis_(upper_bound_axis){};
    virtual ~AlignedBoxShape(){};

    Vecd HalfSize() { return halfsize_; }
    bool checkInBounds(const Vecd &probe_point);
    bool checkUpperBound(const Vecd &probe_point);
    bool checkLowerBound(const Vecd &probe_point);
    bool checkNearUpperBound(const Vecd &probe_point, Real threshold);
    bool checkNearLowerBound(const Vecd &probe_point, Real threshold);
    Vecd getUpperPeriodic(const Vecd &probe_point);
    Vecd getLowerPeriodic(const Vecd &probe_point);
    int AlignmentAxis() { return alignment_axis_; };
};
} // namespace SPH
#endif // COMPLEX_SHAPE_H
