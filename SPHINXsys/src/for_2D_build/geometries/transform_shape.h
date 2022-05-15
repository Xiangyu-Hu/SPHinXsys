/**
 * @file multi_polygon_shape.h
 * @brief Here, we define the 2D geometric algortihms. they are based on the boost library.
 * @details The idea is to define complex geometry based on shapes, usually
 * multi-polygon using boost library. we propose only very simple combinaton
 * that the region is composed of shapes without intersection.
 * That is, the shapes are those contain each other or without overlap.
 * This strict requirement suggests that complex shapes should be finished
 * already in modeling using related binary operations before it is included.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#ifndef TRANSFORM_2d_SHAPE_H
#define TRANSFORM_2d_SHAPE_H

#include "base_data_package.h"
#include "base_geometry.h"

namespace SPH
{
    /**
     * @class TransformShape
     * @brief A template shape in which coordinate transformation is applied applied befone access the interface functions.
     */
    template <class BaseShapeType>
    class TransformShape : public BaseShapeType
    {

    public:
        /** Default constructor. */
        template <typename... ConstructorArgs>
        explicit TransformShape(const Transform2d &transform_2d, ConstructorArgs &&...args)
            : BaseShapeType(std::forward<ConstructorArgs>(args)...),
              transform_2d_(transform_2d), rotation_(transform_2d_.Rotation()){};
        virtual ~TransformShape(){};

        virtual BoundingBox findBounds() override
        {
            BoundingBox original_bound = BaseShapeType::findBounds();
            return BoundingBox(transform_2d_.imposeTransform(original_bound.first),
                               transform_2d_.imposeTransform(original_bound.second));
        };
        virtual bool checkContain(const Vec2d &input_pnt, bool BOUNDARY_INCLUDED = true) override
        {
            Vec2d input_pnt_origin = transform_2d_.imposeInverseTransform(input_pnt);
            return BaseShapeType::checkContain(input_pnt_origin);
        };
        virtual Vec2d findClosestPoint(const Vec2d &input_pnt) override
        {
            Vec2d input_pnt_origin = transform_2d_.imposeInverseTransform(input_pnt);
            Vec2d closest_point_origin = BaseShapeType::findClosestPoint(input_pnt_origin);
            return transform_2d_.imposeTransform(closest_point_origin);
        };
        virtual bool checkNotFar(const Vec2d &input_pnt, Real threshold) override
        {
            Vec2d input_pnt_origin = transform_2d_.imposeInverseTransform(input_pnt);
            return BaseShapeType::checkNotFar(input_pnt_origin, threshold);
        };
        virtual bool checkNearSurface(const Vec2d &input_pnt, Real threshold) override
        {
            Vec2d input_pnt_origin = transform_2d_.imposeInverseTransform(input_pnt);
            return BaseShapeType::checkNearSurface(input_pnt_origin, threshold);
        };
        virtual Real findSignedDistance(const Vec2d &input_pnt) override
        {
            Vec2d input_pnt_origin = transform_2d_.imposeInverseTransform(input_pnt);
            return BaseShapeType::findSignedDistance(input_pnt_origin);
        };
        virtual Vec2d findNormalDirection(const Vec2d &input_pnt) override
        {
            Vec2d input_pnt_origin = transform_2d_.imposeInverseTransform(input_pnt);
            Vec2d normal_direction_origin = BaseShapeType::findNormalDirection(input_pnt_origin);
            return rotation_.imposeTransform(normal_direction_origin);
        };

    protected:
        Transform2d transform_2d_;
        Rotation2d rotation_;
    };
}

#endif // TRANSFORM_2d_SHAPE_H