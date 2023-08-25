#include "geometric_shape.h"

namespace SPH
{
//=================================================================================================//
bool GeometricShape::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
{
    SimTK::UnitVec3 normal;
    bool inside = false;
    contact_geometry_->findNearestPoint(SimTKVec3(probe_point[0], probe_point[1], probe_point[2]), inside, normal);

    return inside;
}
//=================================================================================================//
Vec3d GeometricShape::findClosestPoint(const Vec3d &probe_point)
{
    SimTK::UnitVec3 normal;
    bool inside = false;
    SimTKVec3 out_pnt = contact_geometry_->findNearestPoint(SimTKVec3(probe_point[0], probe_point[1], probe_point[2]), inside, normal);

    return Vecd(out_pnt[0], out_pnt[1], out_pnt[2]);
}
//=================================================================================================//
GeometricShapeBox::
    GeometricShapeBox(const Vecd &halfsize, const std::string &shape_name)
    : GeometricShape(shape_name), brick_(EigenToSimTK(halfsize)), halfsize_(halfsize)
{
    contact_geometry_ = &brick_;
}
//=================================================================================================//
bool GeometricShapeBox::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
{
    return brick_.getGeoBox().containsPoint(SimTKVec3(probe_point[0], probe_point[1], probe_point[2]));
}
//=================================================================================================//
Vec3d GeometricShapeBox::findClosestPoint(const Vec3d &probe_point)
{
    bool inside = false;
    SimTKVec3 out_pnt = brick_.getGeoBox().findClosestPointOnSurface(SimTKVec3(probe_point[0], probe_point[1], probe_point[2]), inside);

    return Vecd(out_pnt[0], out_pnt[1], out_pnt[2]);
}
//=================================================================================================//
BoundingBox GeometricShapeBox::findBounds()
{
    return BoundingBox(-halfsize_, halfsize_);
}
//=================================================================================================//
GeometricShapeBall::
    GeometricShapeBall(const Vecd &center, const Real &radius, const std::string &shape_name)
    : GeometricShape(shape_name), center_(center), sphere_(radius)
{
    contact_geometry_ = &sphere_;
}
//=================================================================================================//
bool GeometricShapeBall::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
{
    return (probe_point - center_).norm() < sphere_.getRadius();
}
//=================================================================================================//
Vec3d GeometricShapeBall::findClosestPoint(const Vec3d &probe_point)
{
    Vec3d displacement = probe_point - center_;
    Real distance = displacement.norm();
    Real cosine0 = (SGN(displacement[0]) * (ABS(displacement[0])) + TinyReal) / (distance + TinyReal);
    Real cosine1 = displacement[1] / (distance + TinyReal);
    Real cosine2 = displacement[2] / (distance + TinyReal);
    return probe_point + (sphere_.getRadius() - distance) * Vec3d(cosine0, cosine1, cosine2);
}
//=================================================================================================//
BoundingBox GeometricShapeBall::findBounds()
{
    Vecd shift = Vecd(sphere_.getRadius(), sphere_.getRadius(), sphere_.getRadius());
    return BoundingBox(center_ - shift, center_ + shift);
}
//=================================================================================================//
} // namespace SPH