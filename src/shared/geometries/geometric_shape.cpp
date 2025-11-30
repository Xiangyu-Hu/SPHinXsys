#include "geometric_shape.h"

namespace SPH
{
//=================================================================================================//
GeometricShapeBox::GeometricShapeBox(
    const Transform &transform, const Vecd &halfsize, const std::string &name)
    : TransformShape<GeometricBox>(name, transform, halfsize) {}
//=================================================================================================//
GeometricShapeBox::GeometricShapeBox(const BoundingBoxd &bounding_box, const std::string &name)
    : TransformShape<GeometricBox>(
          name,
          Transform(0.5 * (bounding_box.lower_ + bounding_box.upper_)),
          0.5 * (bounding_box.upper_ - bounding_box.lower_)) {}
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
BoundingBoxd GeometricShapeBall::findBounds()
{
    Vecd shift = radius_ * Vecd::Ones();
    return BoundingBoxd(center_ - shift, center_ + shift);
}
//=================================================================================================//
} // namespace SPH