#ifndef TRANSFORM_SHAPE_HPP
#define TRANSFORM_SHAPE_HPP

#include "transform_geometry.h"

namespace SPH
{
//=================================================================================================//
template <class GeometryType>
BoundingBoxd TransformGeometry<GeometryType>::findBounds()
{
    BoundingBoxd original_bound = GeometryType::findBounds();
    Vec3d bb_min = Vec3d::Constant(MaxReal);
    Vec3d bb_max = Vec3d::Constant(-MaxReal);
    for (auto x : {original_bound.lower_.x(), original_bound.upper_.x()})
    {
        for (auto y : {original_bound.lower_.y(), original_bound.upper_.y()})
        {
            for (auto z : {original_bound.lower_.z(), original_bound.upper_.z()})
            {
                bb_min = bb_min.cwiseMin(this->transform_.shiftFrameStationToBase(Vec3d(x, y, z)));
                bb_max = bb_max.cwiseMax(this->transform_.shiftFrameStationToBase(Vec3d(x, y, z)));
            }
        }
    }
    return BoundingBoxd(bb_min, bb_max);
}
//=================================================================================================//
} // namespace SPH

#endif // TRANSFORM_SHAPE_HPP