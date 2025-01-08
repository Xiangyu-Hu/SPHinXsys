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
 * @file 	transform_shape.h
 * @brief 	transformation related class for geometries.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef TRANSFORM_SHAPE_H
#define TRANSFORM_SHAPE_H

#include "base_data_package.h"
#include "base_geometry.h"

namespace SPH
{
/**
 * @class TransformShape
 * @brief A template shape in which coordinate transformation is applied
 * before or/and after access the interface functions.
 * Note that this is more suitable to apply for simple geometric shapes.
 */
template <class BaseShapeType>
class TransformShape : public BaseShapeType
{

  public:
    /** template constructor for general shapes. */
    template <typename... Args>
    explicit TransformShape(const Transform &transform, Args &&...args)
        : BaseShapeType(std::forward<Args>(args)...), transform_(transform){};

    virtual ~TransformShape(){};

    /** variable transform is introduced here */
    Transform &getTransform() { return transform_; };
    void setTransform(const Transform &transform) { transform_ = transform; };

    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override
    {
        Vecd input_pnt_origin = transform_.shiftBaseStationToFrame(probe_point);
        return BaseShapeType::checkContain(input_pnt_origin);
    };

    virtual Vecd findClosestPoint(const Vecd &probe_point) override
    {
        Vecd input_pnt_origin = transform_.shiftBaseStationToFrame(probe_point);
        Vecd closest_point_origin = BaseShapeType::findClosestPoint(input_pnt_origin);
        return transform_.shiftFrameStationToBase(closest_point_origin);
    };

  protected:
    Transform transform_;

    /// Returns the AABB of the rotated underlying shape's AABB
    /// ⚠️ It is not the tight fit AABB of the underlying shape
    /// But at least it encloses the underlying shape fully
    virtual BoundingBox findBounds() override
    {
        BoundingBox original_bound = BaseShapeType::findBounds();
        Vecd bb_min = Vecd::Constant(MaxReal);
        Vecd bb_max = Vecd::Constant(-MaxReal);
        for (auto x : {original_bound.first_.x(), original_bound.second_.x()})
        {
            for (auto y : {original_bound.first_.y(), original_bound.second_.y()})
            {
                if constexpr (Dimensions == 3)
                {
                    for (auto z : {original_bound.first_.z(), original_bound.second_.z()})
                    {
                        bb_min = bb_min.cwiseMin(transform_.shiftFrameStationToBase(Vecd(x, y, z)));
                        bb_max = bb_max.cwiseMax(transform_.shiftFrameStationToBase(Vecd(x, y, z)));
                    }
                }
                else
                {
                    bb_min = bb_min.cwiseMin(transform_.shiftFrameStationToBase(Vecd(x, y)));
                    bb_max = bb_max.cwiseMax(transform_.shiftFrameStationToBase(Vecd(x, y)));
                }
            }
        }
        return BoundingBox(bb_min, bb_max);
    };
};
} // namespace SPH

#endif // TRANSFORM_SHAPE_H
