/**
 * @file transform_shape.h
 * @brief tranformation related class for geometries.
 * @author	Xiangyu Hu
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
     */
    template <class BaseShapeType>
    class TransformShape : public BaseShapeType
    {

    public:
        /** template constructor for general shapes. */
        template <typename... ConstructorArgs>
        explicit TransformShape(const Transformd &transformd, ConstructorArgs &&...args)
            : BaseShapeType(std::forward<ConstructorArgs>(args)...),
              transformd_(transformd){};
        virtual ~TransformShape(){};

        /** variable transform is introduced here */
        Transformd getTransform() { return transformd_; };
        void setTransform(const Transformd &transformd) { transformd_ = transformd; };

        virtual BoundingBox findBounds() override
        {
            BoundingBox original_bound = BaseShapeType::findBounds();
            return BoundingBox(transformd_.shiftFrameStationToBase(original_bound.first),
                               transformd_.shiftFrameStationToBase(original_bound.second));
        };
        virtual bool checkContain(const Vecd &input_pnt, bool BOUNDARY_INCLUDED = true) override
        {
            Vecd input_pnt_origin = transformd_.shiftBaseStationToFrame(input_pnt);
            return BaseShapeType::checkContain(input_pnt_origin);
        };
        virtual Vecd findClosestPoint(const Vecd &input_pnt) override
        {
            Vecd input_pnt_origin = transformd_.shiftBaseStationToFrame(input_pnt);
            Vecd closest_point_origin = BaseShapeType::findClosestPoint(input_pnt_origin);
            return transformd_.shiftFrameStationToBase(closest_point_origin);
        };
        virtual bool checkNotFar(const Vecd &input_pnt, Real threshold) override
        {
            Vecd input_pnt_origin = transformd_.shiftBaseStationToFrame(input_pnt);
            return BaseShapeType::checkNotFar(input_pnt_origin, threshold);
        };
        virtual bool checkNearSurface(const Vecd &input_pnt, Real threshold) override
        {
            Vecd input_pnt_origin = transformd_.shiftBaseStationToFrame(input_pnt);
            return BaseShapeType::checkNearSurface(input_pnt_origin, threshold);
        };
        virtual Real findSignedDistance(const Vecd &input_pnt) override
        {
            Vecd input_pnt_origin = transformd_.shiftBaseStationToFrame(input_pnt);
            return BaseShapeType::findSignedDistance(input_pnt_origin);
        };
        virtual Vecd findNormalDirection(const Vecd &input_pnt) override
        {
            Vecd input_pnt_origin = transformd_.shiftBaseStationToFrame(input_pnt);
            Vecd normal_direction_origin = BaseShapeType::findNormalDirection(input_pnt_origin);
            return transformd_.xformFrameVecToBase(normal_direction_origin);
        };

    protected:
        Transformd transformd_;
    };
}

#endif // TRANSFORM_SHAPE_H