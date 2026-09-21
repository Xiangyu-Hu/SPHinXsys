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

#include "base_geometry.h"
#include "data_type.h"

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
        : BaseShapeType(std::forward<Args>(args)...) {};

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
class ExtrudeShape : public Shape
{
    SharedPtrKeeper<Shape> base_shape_keeper_;
    Real thickness_, thickness_sqr_;
    Shape &base_shape_;
    Vecd getShift(const Vecd &probe_point, const Vecd &original_closest_point);

  public:
    template <class BaseShapeType>
    ExtrudeShape(Real thickness, SharedPtr<BaseShapeType> base_shape_ptr)
        : Shape("Extruded" + base_shape_ptr->Name()),
          thickness_(thickness), thickness_sqr_(thickness * thickness),
          base_shape_(*base_shape_keeper_.assignPtr(base_shape_ptr)) {};
    ExtrudeShape(Shape &base_shape, Real thickness, const std::string &shape_name = "ExtrudedShape");
    virtual ~ExtrudeShape() {};
    virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override;
    virtual Vecd findClosestPoint(const Vecd &probe_point) override;
    virtual BoundingBoxd findBounds() override;
};
} // namespace SPH

#endif // MAPPING_SHAPE_H
