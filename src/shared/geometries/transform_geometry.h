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
 * @file 	transform_geometry.h
 * @brief 	Linear transformation related class for geometries.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef TRANSFORM_SHAPE_H
#define TRANSFORM_SHAPE_H

#include "base_data_type_package.h"
#include "base_geometry.h"

namespace SPH
{
template <class GeometryType>
class TransformGeometry : public GeometryType
{

  public:
    template <typename... Args>
    explicit TransformGeometry(const Transform &transform, Args &&...args)
        : GeometryType(std::forward<Args>(args)...), transform_(transform){};
    ~TransformGeometry() {};

    /** variable transform is introduced here */
    Transform &getTransform() { return transform_; };
    void setTransform(const Transform &transform) { transform_ = transform; };

    bool checkContain(const Vecd &probe_point)
    {
        Vecd input_pnt_origin = transform_.shiftBaseStationToFrame(probe_point);
        return GeometryType::checkContain(input_pnt_origin);
    };

    Vecd findClosestPoint(const Vecd &probe_point)
    {
        Vecd input_pnt_origin = transform_.shiftBaseStationToFrame(probe_point);
        Vecd closest_point_origin = GeometryType::findClosestPoint(input_pnt_origin);
        return transform_.shiftFrameStationToBase(closest_point_origin);
    };

  protected:
    Transform transform_;
};

/**
 * @class TransformShape
 * @brief A template shape in which coordinate transformation is applied
 * before or/and after access the interface functions.
 * Note that this is more suitable to apply for simple geometric shapes.
 */
template <class GeometryType>
class TransformShape : public TransformGeometry<GeometryType>, public Shape
{

  public:
    /** template constructor for general shapes. */
    template <typename... Args>
    explicit TransformShape(const std::string &name, const Transform &transform, Args &&...args)
        : TransformGeometry<GeometryType>(transform, std::forward<Args>(args)...),
          Shape(name){};
    virtual ~TransformShape() {};

    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override
    {
        return TransformGeometry<GeometryType>::checkContain(probe_point);
    };

    virtual Vecd findClosestPoint(const Vecd &probe_point) override
    {
        return TransformGeometry<GeometryType>::findClosestPoint(probe_point);
    };

    // Returns the AABB of the rotated underlying shape's AABB
    // It is not the tight fit AABB of the underlying shape
    // But at least it encloses the underlying shape fully
    virtual BoundingBoxd findBounds() override
    {
        BoundingBoxd original_bound = TransformGeometry<GeometryType>::findBounds();
        Vecd bb_min = Vecd::Constant(MaxReal);
        Vecd bb_max = Vecd::Constant(-MaxReal);
        for (auto x : {original_bound.lower_.x(), original_bound.upper_.x()})
        {
            for (auto y : {original_bound.lower_.y(), original_bound.upper_.y()})
            {
                if constexpr (Dimensions == 3)
                {
                    for (auto z : {original_bound.lower_.z(), original_bound.upper_.z()})
                    {
                        bb_min = bb_min.cwiseMin(this->transform_.shiftFrameStationToBase(Vecd(x, y, z)));
                        bb_max = bb_max.cwiseMax(this->transform_.shiftFrameStationToBase(Vecd(x, y, z)));
                    }
                }
                else
                {
                    bb_min = bb_min.cwiseMin(this->transform_.shiftFrameStationToBase(Vecd(x, y)));
                    bb_max = bb_max.cwiseMax(this->transform_.shiftFrameStationToBase(Vecd(x, y)));
                }
            }
        }
        return BoundingBoxd(bb_min, bb_max);
    };
};
} // namespace SPH

#endif // TRANSFORM_SHAPE_H
