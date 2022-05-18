#include "geometric_shape.h"

namespace SPH
{
    //=================================================================================================//
    bool GeometricShape::checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        contact_geometry_->findNearestPoint(pnt, inside, normal);

        return inside;
    }
    //=================================================================================================//
    Vec3d GeometricShape::findClosestPoint(const Vec3d &pnt)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        return contact_geometry_->findNearestPoint(pnt, inside, normal);
    }
    //=================================================================================================//
    GeometricShapeBrick::
        GeometricShapeBrick(const Vec3d &halfsize, const std::string &shape_name)
        : GeometricShape(shape_name), brick_(halfsize), halfsize_(halfsize)
    {
        contact_geometry_ = &brick_;
    }
    //=================================================================================================//
    bool GeometricShapeBrick::checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED)
    {
        return brick_.getGeoBox().containsPoint(pnt);
    }
    //=================================================================================================//
    Vec3d GeometricShapeBrick::findClosestPoint(const Vec3d &pnt)
    {
        bool inside = false;
        return brick_.getGeoBox().findClosestPointOnSurface(pnt, inside);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeBrick::findBounds()
    {
         return BoundingBox(- halfsize_, halfsize_);
    }
    //=================================================================================================//
    GeometricShapeSphere::
        GeometricShapeSphere(const Vec3d &center, const Real &radius, const std::string &shape_name)
        : GeometricShape(shape_name), center_(center), sphere_(radius)
    {
        contact_geometry_ = &sphere_;
    }
    //=================================================================================================//
    bool GeometricShapeSphere::checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED)
    {
        return (pnt - center_).norm() < sphere_.getRadius();
    }
    //=================================================================================================//
    Vec3d GeometricShapeSphere::findClosestPoint(const Vec3d &pnt)
    {
        Vec3d displacement = pnt - center_;
        Real distance = displacement.norm();
        Real cosine0 = (SGN(displacement[0]) * (ABS(displacement[0])) + TinyReal) / (distance + TinyReal);
        Real cosine1 = displacement[1] / (distance + TinyReal);
        Real cosine2 = displacement[2] / (distance + TinyReal);
        return pnt + (sphere_.getRadius() - distance) * Vec3d(cosine0, cosine1, cosine2);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeSphere::findBounds()
    {
        Vec3d shift = Vec3d(sphere_.getRadius(), sphere_.getRadius(), sphere_.getRadius());
        return BoundingBox(center_ - shift, center_ + shift);
    }
    //=================================================================================================//
}