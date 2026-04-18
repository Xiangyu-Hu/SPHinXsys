#ifndef SDF_OPERATION_HPP
#define SDF_OPERATION_HPP

#include "sdf_operation.h"

namespace SPH
{
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFAddition::operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    return SMIN(input1(point), input2(point));
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFAddition::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().add(input2.findBounds());
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSubtraction::operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    return SMAX(input1(point), -input2(point));
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFIntersection::operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    return SMAX(input1(point), input2(point));
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFIntersection::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().intersect(input2.findBounds());
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSmoothAddition::operator()(
    const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    Real d1 = input1(point);
    Real d2 = input2(point);
    Real k = 4.0 * finest_grid_spacing_;
    Real h = SMAX(k - ABS(d1 - d2), Real(0));
    return SMIN(d1, d2) - h * h * 0.25 / k;
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFSmoothAddition::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().add(input2.findBounds());
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSmoothSubtraction::operator()(
    const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    Real d1 = -input1(point);
    Real d2 = input2(point);
    Real k = 4.0 * finest_grid_spacing_;
    Real h = SMAX(k - ABS(d1 - d2), Real(0));
    return -(SMIN(d1, d2) - h * h * 0.25 / k);
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFSmoothSubtraction::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds();
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSmoothIntersection::operator()(
    const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    Real d1 = -input1(point);
    Real d2 = -input2(point);
    Real k = 4.0 * finest_grid_spacing_;
    Real h = SMAX(k - ABS(d1 - d2), Real(0));
    return -(SMIN(d1, d2) - h * h * 0.25 / k);
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFSmoothIntersection::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().intersect(input2.findBounds());
}
//=================================================================================================//
template <typename Input2D>
Real SDFExtrusion::operator()(const Input2D &input, const Vec3d &point) const
{
    Vec2d point_2d = point.tail(2); // Project to 2D plane
    Real axial = point[0];
    if (axial < 0.0 || axial > height_)
        return MaxReal; // Outside the extrusion's height
    return input(point_2d);
}
//=================================================================================================//
template <typename Input2D>
Real SDFRotation::operator()(const Input2D &input, const Vec3d &point) const
{
    Real cos_angle = std::cos(angle_);
    Real sin_angle = std::sin(angle_);
    Vec2d rotated_point;
    rotated_point[0] = cos_angle * point[0] - sin_angle * point[1];
    rotated_point[1] = sin_angle * point[0] + cos_angle * point[1];
    return input(rotated_point);
}
//=================================================================================================//
template <typename Input3D>
Real SDFElongation::operator()(const Input3D &input, const Vec3d &point) const
{
    Vec3d elongated_point = point;
    elongated_point[0] *= elongation_factor_; // Elongate along x-axis
    return input(elongated_point);
}
//=================================================================================================//
} // namespace SPH
#endif // SDF_OPERATION_HPP