#include "geometric_shape.h"

namespace SPH
{
    //=================================================================================================//
    bool GeometricShape::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        contact_geometry_->findNearestPoint(SimTK::Vec3(pnt[0], pnt[1], pnt[2]), inside, normal);

        return inside;
    }
    //=================================================================================================//
    Vecd GeometricShape::findClosestPoint(const Vecd &pnt)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        SimTK::Vec3 out_pnt = contact_geometry_->findNearestPoint(SimTK::Vec3(pnt[0], pnt[1], pnt[2]), inside, normal);

        return Vecd(out_pnt[0], out_pnt[1], out_pnt[2]);
    }
    //=================================================================================================//
    GeometricShapeBox::
        GeometricShapeBox(const Vecd &halfsize, const std::string &shape_name)
        : GeometricShape(shape_name), brick_(halfsize), halfsize_(halfsize)
    {
        contact_geometry_ = &brick_;
    }
    //=================================================================================================//
    bool GeometricShapeBox::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED)
    {
        return brick_.getGeoBox().containsPoint(SimTK::Vec3(pnt[0], pnt[1], pnt[2]));
    }
    //=================================================================================================//
    Vecd GeometricShapeBox::findClosestPoint(const Vecd &pnt)
    {
        bool inside = false;
        SimTK::Vec3 out_pnt = brick_.getGeoBox().findClosestPointOnSurface(SimTK::Vec3(pnt[0], pnt[1], pnt[2]), inside);

        return Vecd(out_pnt[0], out_pnt[1], out_pnt[2]);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeBox::findBounds()
    {
         return BoundingBox(- halfsize_, halfsize_);
    }
    //=================================================================================================//
    GeometricShapeBall::
        GeometricShapeBall(const Vecd &center, const Real &radius, const std::string &shape_name)
        : GeometricShape(shape_name), center_(center), sphere_(radius)
    {
        contact_geometry_ = &sphere_;
    }
    //=================================================================================================//
    bool GeometricShapeBall::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED)
    {
        return (pnt - center_).norm() < sphere_.getRadius();
    }
    //=================================================================================================//
    Vecd GeometricShapeBall::findClosestPoint(const Vecd &pnt)
    {
        Vecd displacement = pnt - center_;
        Real distance = displacement.norm();
        Real cosine0 = (SGN(displacement[0]) * (ABS(displacement[0])) + TinyReal) / (distance + TinyReal);
        Real cosine1 = displacement[1] / (distance + TinyReal);
        Real cosine2 = displacement[2] / (distance + TinyReal);
        return pnt + (sphere_.getRadius() - distance) * Vecd(cosine0, cosine1, cosine2);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeBall::findBounds()
    {
        Vecd shift = Vecd(sphere_.getRadius(), sphere_.getRadius(), sphere_.getRadius());
        return BoundingBox(center_ - shift, center_ + shift);
    }
    //=================================================================================================//
}