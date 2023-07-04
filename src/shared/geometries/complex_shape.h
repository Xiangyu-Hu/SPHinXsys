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

    template <typename... ConstructorArgs>
    LevelSetShape *defineLevelSetShape(SPHBody &sph_body, const std::string &shape_name, ConstructorArgs &&...args)
    {
        size_t index = getShapeIndexByName(shape_name);
        LevelSetShape *level_set_shape = shapes_ptr_keeper_[index].createPtr<LevelSetShape>(
            sph_body, *shapes_and_ops_[index].first, std::forward<ConstructorArgs>(args)...);
        shapes_and_ops_[index].first = DynamicCast<Shape>(this, level_set_shape);
        return level_set_shape;
    };
};

using DefaultShape = ComplexShape;

/**
 * @class AlignedBoxShape
 * @brief Used to describe a bounding box in which
 * the plane vertical to axis direction is aligned to a planar piece of a shape.
 */
class AlignedBoxShape : public TransformShape<GeometricShapeBox>
{
  public:
    /** construct directly */
    template <typename... Args>
    explicit AlignedBoxShape(const Transform &transform, Args &&...args)
        : TransformShape<GeometricShapeBox>(transform, std::forward<Args>(args)...){};
    /** construct from a shape already has aligned boundaries */
    template <typename... Args>
    explicit AlignedBoxShape(const Shape &shape, Args &&...args)
        : TransformShape<GeometricShapeBox>(
              Transform(Vecd(0.5 * (shape.bounding_box_.second_ + shape.bounding_box_.first_))),
              0.5 * (shape.bounding_box_.second_ - shape.bounding_box_.first_), std::forward<Args>(args)...){};

    Vecd HalfSize() { return halfsize_; }
    bool checkInBounds(int axis, const Vecd &probe_point);
    bool checkUpperBound(int axis, const Vecd &probe_point);
    bool checkLowerBound(int axis, const Vecd &probe_point);
    bool checkNearUpperBound(int axis, const Vecd &probe_point, Real threshold);
    bool checkNearLowerBound(int axis, const Vecd &probe_point, Real threshold);
    Vecd getUpperPeriodic(int axis, const Vecd &probe_point);
    Vecd getLowerPeriodic(int axis, const Vecd &probe_point);
};
} // namespace SPH

#endif // COMPLEX_SHAPE_H