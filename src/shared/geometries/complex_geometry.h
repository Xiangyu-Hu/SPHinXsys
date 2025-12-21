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
 * @file complex_geometry.h
 * @brief Here, we define the a container of different type of shapes, which allow
 * the boolean operation between the different type of shapes. The shapes
 * can be defined previously and add to this complex shapes container.
 */

#ifndef COMPLEX_SHAPE_H
#define COMPLEX_SHAPE_H

#include "base_geometry.h"
#include "geometric_shape.h"
#include "level_set_shape.h"
#include "transform_geometry.h"

namespace SPH
{
class LevelSetShape;

/**
 * @class ComplexShape
 * @brief  For now, if the level set shape (for particle relaxation)
 * will be generated from the complex shape,
 * partially overlapped shapes are not allowed for 3D problems.
 * However, if only the contain function
 * is used, for example generating particles using lattice generator,
 * partially overlapped shapes are allowed.
 **/
class ComplexShape : public BinaryShapes
{
  public:
    explicit ComplexShape(const std::string &shape_name)
        : BinaryShapes(shape_name) {};
    virtual ~ComplexShape() {};

    template <typename... Args>
    LevelSetShape *defineLevelSetShape(SPHBody &sph_body, const std::string &shape_name, Args &&...args)
    {
        size_t index = getSubShapeIndexByName(shape_name);
        LevelSetShape *level_set_shape = sub_shapes_keeper_[index].createPtr<LevelSetShape>(
            sph_body, *sub_shapes_and_ops_[index].first, std::forward<Args>(args)...);
        sub_shapes_and_ops_[index].first = DynamicCast<Shape>(this, level_set_shape);
        return level_set_shape;
    };
};

using DefaultShape = ComplexShape;

/**
 * @class AlignedBox
 * @brief Used to describe a bounding box in which
 * the upper bound direction is aligned to the normal of a planar piece on a shape.
 */
class AlignedBox : public TransformGeometry<GeometricBox>
{
    int alignment_axis_;

  public:
    /** construct directly */
    template <typename... Args>
    explicit AlignedBox(int upper_bound_axis, const Transform &transform, Args &&...args)
        : TransformGeometry<GeometricBox>(transform, std::forward<Args>(args)...),
          alignment_axis_(upper_bound_axis){};
    /** construct from a shape already has aligned boundaries */
    template <typename... Args>
    explicit AlignedBox(int upper_bound_axis, const Shape &shape, Args &&...args)
        : TransformGeometry<GeometricBox>(
              Transform(Vecd(0.5 * (shape.bounding_box_.upper_ + shape.bounding_box_.lower_))),
              0.5 * (shape.bounding_box_.upper_ - shape.bounding_box_.lower_), std::forward<Args>(args)...),
          alignment_axis_(upper_bound_axis){};
    ~AlignedBox() {};

    Vecd HalfSize() { return halfsize_; }
    bool checkNearSurface(const Vecd &probe_point, Real threshold);
    bool checkNotFar(const Vecd &probe_point, Real threshold);
    bool checkInBounds(const Vecd &probe_point, Real lower_bound_fringe = 0.0, Real upper_bound_fringe = 0.0)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return position_in_frame[alignment_axis_] >= -halfsize_[alignment_axis_] - lower_bound_fringe &&
                       position_in_frame[alignment_axis_] <= halfsize_[alignment_axis_] + upper_bound_fringe
                   ? true
                   : false;
    };
    bool checkUpperBound(const Vecd &probe_point, Real upper_bound_fringe = 0.0)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return position_in_frame[alignment_axis_] > halfsize_[alignment_axis_] + upper_bound_fringe ? true : false;
    };
    bool checkLowerBound(const Vecd &probe_point, Real lower_bound_fringe = 0.0)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        return position_in_frame[alignment_axis_] < -halfsize_[alignment_axis_] - lower_bound_fringe ? true : false;
    }
    bool checkNearUpperBound(const Vecd &probe_point, Real threshold);
    bool checkNearLowerBound(const Vecd &probe_point, Real threshold);
    Vecd getUpperPeriodic(const Vecd &probe_point)
    {
        Vecd position_in_frame = transform_.shiftBaseStationToFrame(probe_point);
        Vecd shift = Vecd::Zero();
        shift[alignment_axis_] -= 2.0 * halfsize_[alignment_axis_];
        return transform_.shiftFrameStationToBase(position_in_frame + shift);
    };
    Vecd getLowerPeriodic(const Vecd &probe_point);
    int AlignmentAxis() { return alignment_axis_; };

    void imposeAlignment(Vecd &value)
    {
        Vecd frame_velocity = Vecd::Zero();
        frame_velocity[alignment_axis_] = transform_.xformBaseVecToFrame(value)[alignment_axis_];
        value = transform_.xformFrameVecToBase(frame_velocity);
    };
};
} // namespace SPH
#endif // COMPLEX_SHAPE_H
