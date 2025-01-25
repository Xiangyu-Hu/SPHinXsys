#include "geometric_shape.h"

namespace SPH
{
//=================================================================================================//
GeometricShapeBox::GeometricShapeBox(const Vecd &halfsize, const std::string &shape_name)
    : GeometricBox(halfsize), Shape(shape_name) {}
//=================================================================================================//
bool GeometricShapeBox::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    return GeometricBox::checkContain(probe_point);
}
//=================================================================================================//
Vecd GeometricShapeBox::findClosestPoint(const Vecd &probe_point)
{
    return GeometricBox::findClosestPoint(probe_point);
}
//=================================================================================================//
BoundingBox GeometricShapeBox::findBounds()
{
    return BoundingBox(-halfsize_, halfsize_);
}
//=================================================================================================//
GeometricShapeBall::GeometricShapeBall(const Vecd &center, Real radius,
                                       const std::string &shape_name)
    : GeometricBall(radius), Shape(shape_name), center_(center) {}
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