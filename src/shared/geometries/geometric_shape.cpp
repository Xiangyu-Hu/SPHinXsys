#include "geometric_shape.h"

namespace SPH
{
//=================================================================================================//
GeometricShapeBox::GeometricShapeBox(
    const Transform &transform, const Vecd &halfsize, const std::string &name)
    : TransformShape<GeometricBox>(name, transform, halfsize) {}
//=================================================================================================//
GeometricShapeBox::GeometricShapeBox(const BoundingBox &bounding_box, const std::string &name)
    : TransformShape<GeometricBox>(
          name,
          Transform(0.5 * (bounding_box.first_ + bounding_box.second_)),
          0.5 * (bounding_box.second_ - bounding_box.first_)) {}
//=================================================================================================//
GeometricShapeBox::GeometricShapeBox(
    const Vecd &lower_bound, const Vecd &upper_bound, const std::string &name)
    : TransformShape<GeometricBox>(
          name,
          Transform(0.5 * (lower_bound + upper_bound)),
          0.5 * (upper_bound - lower_bound)) {}
//=================================================================================================//
GeometricShapeBall::GeometricShapeBall(const Vecd &center, Real radius,
                                       const std::string &name)
    : GeometricBall(radius), Shape(name), center_(center) {}
//=================================================================================================//
bool GeometricShapeBall::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    return GeometricBall::checkContain(probe_point - center_);
}
//=================================================================================================//
Vecd GeometricShapeBall::findClosestPoint(const Vecd &probe_point)
{
    return center_ + GeometricBall::findClosestPoint(probe_point - center_);
}
//=================================================================================================//
BoundingBox GeometricShapeBall::findBounds()
{
    Vecd shift = radius_ * Vecd::Ones();
    return BoundingBox(center_ - shift, center_ + shift);
}
//=================================================================================================//
} // namespace SPH