/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	exclusive_shape.h
 * @brief 	exclusive shape related class for geometries.
 * @author	Yongchuan Yu and Xiangyu Hu
 */

#ifndef EXCLUSIVE_SHAPE_H
#define EXCLUSIVE_SHAPE_H

#include "base_data_package.h"
#include "base_geometry.h"

namespace SPH
{
/**
 * @class ExclusiveShape
 * @brief A template shape which has the fluid outside of the geometry.
 */
    class Shape;

    template <class BaseShapeType>
    class ExclusiveShape : public BaseShapeType
    {

      public:
        /** template constructor for general shapes. */
        template <typename... ConstructorArgs>
        explicit ExclusiveShape(ConstructorArgs &&...args)
            : BaseShapeType(std::forward<ConstructorArgs>(args)...), 
            base_shape_(base_shape_keeper_.createPtr<BaseShapeType>(std::forward<ConstructorArgs>(args)...))
        {};

        virtual ~ExclusiveShape(){};

        /*reverse the value of checkContain function*/
        virtual bool checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED = true) override
        {
            bool is_inside = base_shape_->checkContain(probe_point);
            return !is_inside;
        };

        virtual Vecd findClosestPoint(const Vecd &probe_point) override
        {
            return base_shape_->findClosestPoint(probe_point);
        };

      private:
            UniquePtrKeeper<BaseShapeType> base_shape_keeper_;
            Shape *base_shape_;
    };
} // namespace SPH

#endif // TRANSFORM_SHAPE_H