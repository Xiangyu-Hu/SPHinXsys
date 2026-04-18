#ifndef SDF_EXTENSION_HPP
#define SDF_EXTENSION_HPP

#include "sdf_extension.h"

namespace SPH
{
//=================================================================================================//
template <typename InputType, typename ExtensionType>
SDFExtension<InputType, ExtensionType>::SDFExtension(const InputType &input, const ExtensionType &extension)
    : input_(input), extension_(extension) {}
//=================================================================================================//
template <typename InputType>
Real SDFScale::operator()(const InputType &input, const Vec3d &point) const
{
    return input(point / scale_factor_) * scale_factor_;
}
//=================================================================================================//
template <typename InputType>
auto SDFScale::findBounds(const InputType &input) const
{
    return input.findBounds().scale(scale_factor_);
}
//=================================================================================================//
template <typename... Args>
void SDFTransform::setParameters(Args &&...args)
{
    transform_ = Transform3d(std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename InputType>
Real SDFTransform::operator()(const InputType &input, const Vec3d &point) const
{
    Vec3d transformed_point = transform_.shiftBaseStationToFrame(point);
    return input(transformed_point);
}
//=================================================================================================//
template <typename InputType>
auto SDFTransform::findBounds(const InputType &input) const
{
    BoundingBox3d original_bound = input.findBounds();
    Vec3d bb_min = Vec3d::Constant(MaxReal);
    Vec3d bb_max = Vec3d::Constant(-MaxReal);
    for (auto x : {original_bound.lower_.x(), original_bound.upper_.x()})
    {
        for (auto y : {original_bound.lower_.y(), original_bound.upper_.y()})
        {
            for (auto z : {original_bound.lower_.z(), original_bound.upper_.z()})
            {
                bb_min = bb_min.cwiseMin(transform_.shiftFrameStationToBase(Vec3d(x, y, z)));
                bb_max = bb_max.cwiseMax(transform_.shiftFrameStationToBase(Vec3d(x, y, z)));
            }
        }
    }
    return BoundingBox3d(bb_min, bb_max);
}
//=================================================================================================//
} // namespace SPH
#endif // SDF_EXTENSION_HPP