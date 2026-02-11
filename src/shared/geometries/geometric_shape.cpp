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
GeometricShapeCylinder::GeometricShapeCylinder(const Vecd &center, const Vecd &axis, Real radius, Real halflength,
                                               const std::string &name)
    : GeometricCylinder(axis, radius, halflength), Shape(name), center_(center) {}
//=================================================================================================//
bool GeometricShapeCylinder::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    return GeometricCylinder::checkContain(probe_point - center_);
}
//=================================================================================================//
Vecd GeometricShapeCylinder::findClosestPoint(const Vecd &probe_point)
{
    return center_ + GeometricCylinder::findClosestPoint(probe_point - center_);
}
//=================================================================================================//
BoundingBoxd GeometricShapeCylinder::findBounds()
{
    // Get bounds from base class (centered at origin)
    BoundingBoxd base_bounds = GeometricCylinder::findBounds();
    // Translate to actual center
    return BoundingBoxd(center_ + base_bounds.lower_, center_ + base_bounds.upper_);
}
//=================================================================================================//
} // namespace SPH