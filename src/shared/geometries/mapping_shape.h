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
 * @file 	maping_shape.h
 * @brief 	Shape generated based on mapping from an original shape.
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#ifndef MAPPING_SHAPE_H
#define MAPPING_SHAPE_H

#include "base_data_package.h"
#include "base_geometry.h"

namespace SPH
{
/**
 * @class InverseShape
 * @brief A template shape which switches the defined inside or outside region.
 * @brief In simple terms, it gives opposite return value for the function checkContain() as the original shape
 */

template <class BaseShapeType>
class InverseShape : public BaseShapeType
{

  public:
    /** template constructor for general shapes. */
    template <typename... Args>
    explicit InverseShape(Args &&...args)
        : BaseShapeType(std::forward<Args>(args)...){};

    virtual ~InverseShape() {};

    /*reverse the value of checkContain function*/
    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override
    {
        return !BaseShapeType::checkContain(probe_point);
    };
};

/**
 * @class ExtrudeShape
 * @brief A template shape which define the region by expanding the geometry surface with given thickness.
 * @brief Positive thickness will extend the shape and negative thickness will shrink the shape.
 */
template <class BaseShapeType>
class ExtrudeShape : public Shape
{
    UniquePtrKeeper<BaseShapeType> base_shape_keeper_;
    Real thickness_, thickness_sqr_;

    BaseShapeType *base_shape_;
    Vecd getShift(const Vecd &probe_point, const Vecd &original_closest_point)
    {
        Vecd displacement = original_closest_point - probe_point;
        return thickness_ * displacement / (displacement.norm() + Eps);
    };

  public:
    explicit ExtrudeShape(BaseShapeType *base_shape, Real thickness)
        : Shape("Extruded" + base_shape->getName()),
          thickness_(thickness), thickness_sqr_(thickness * thickness),
          base_shape_(base_shape) {};

    template <typename... Args>
    explicit ExtrudeShape(Real thickness, Args &&...args)
        : Shape("Extruded"), thickness_(thickness), thickness_sqr_(thickness * thickness),
          base_shape_(base_shape_keeper_.template createPtr<BaseShapeType>(std::forward<Args>(args)...))
    {
        name_ = "Extruded" + base_shape_->getName();
    };
    virtual ~ExtrudeShape() {};

    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override
    {
        Vecd original_closest_point = base_shape_->findClosestPoint(probe_point);
        Vecd displacement = original_closest_point - probe_point;
        if (base_shape_->checkContain(probe_point))
        {
            return thickness_ > 0.0 ? true : displacement.squaredNorm() > thickness_sqr_;
        }
        else
        {
            return thickness_ < 0.0 ? false : displacement.squaredNorm() < thickness_sqr_;
        }
    };

    virtual Vecd findClosestPoint(const Vecd &probe_point) override
    {
        Vecd closest_point = base_shape_->findClosestPoint(probe_point);
        Vecd shift = getShift(probe_point, closest_point);
        closest_point += base_shape_->checkContain(probe_point) ? shift : -shift;
        return closest_point;
    };

    virtual BoundingBox findBounds() override
    {
        BoundingBox bounds = base_shape_->findBounds();
        return bounds.expand(thickness_ * Vecd::Ones());
    };
};
} // namespace SPH

#endif // MAPPING_SHAPE_H
