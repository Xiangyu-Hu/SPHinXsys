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
    Vec2d bb_min = Vec2d::Constant(MaxReal);
    Vec2d bb_max = Vec2d::Constant(-MaxReal);
    for (auto x : {original_bound.lower_.x(), original_bound.upper_.x()})
    {
        for (auto y : {original_bound.lower_.y(), original_bound.upper_.y()})
        {
            bb_min = bb_min.cwiseMin(this->transform_.shiftFrameStationToBase(Vec2d(x, y)));
            bb_max = bb_max.cwiseMax(this->transform_.shiftFrameStationToBase(Vec2d(x, y)));
        }
    }
    return BoundingBoxd(bb_min, bb_max);
}
//=================================================================================================//

} // namespace SPH

#endif // TRANSFORM_SHAPE_HPP