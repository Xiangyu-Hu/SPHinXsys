#include "geometric_shape.h"

namespace SPH
{
    //=================================================================================================//
    bool GeometricShape::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        contact_geometry_->findNearestPoint(probe_point, inside, normal);

        return inside;
    }
    //=================================================================================================//
    Vec3d GeometricShape::findClosestPoint(const Vec3d &probe_point)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        return contact_geometry_->findNearestPoint(probe_point, inside, normal);
    }
    //=================================================================================================//
    GeometricShapeBox::
        GeometricShapeBox(const Vec3d &halfsize, const std::string &shape_name)
        : GeometricShape(shape_name), brick_(halfsize), halfsize_(halfsize)
    {
        contact_geometry_ = &brick_;
    }
    //=================================================================================================//
    bool GeometricShapeBox::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
    {
        return brick_.getGeoBox().containsPoint(probe_point);
    }
    //=================================================================================================//
    Vec3d GeometricShapeBox::findClosestPoint(const Vec3d &probe_point)
    {
        bool inside = false;
        return brick_.getGeoBox().findClosestPointOnSurface(probe_point, inside);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeBox::findBounds()
    {
         return BoundingBox(- halfsize_, halfsize_);
    }
    //=================================================================================================//
    GeometricShapeBall::
        GeometricShapeBall(const Vec3d &center, const Real &radius, const std::string &shape_name)
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
        Vec3d shift = Vec3d(sphere_.getRadius(), sphere_.getRadius(), sphere_.getRadius());
        return BoundingBox(center_ - shift, center_ + shift);
    }
    //=================================================================================================//
}