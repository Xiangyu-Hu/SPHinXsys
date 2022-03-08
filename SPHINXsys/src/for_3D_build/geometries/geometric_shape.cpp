#include "geometric_shape.h"

namespace SPH
{
    //=================================================================================================//
    bool GeometricShape::checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        contact_geometry_->findNearestPoint(transform_.shiftBaseStationToFrame(pnt), inside, normal);

        return inside;
    }
    //=================================================================================================//
    Vec3d GeometricShape::findClosestPoint(const Vec3d &pnt)
    {
        SimTK::UnitVec3 normal;
        bool inside = false;
        Vec3d closest_point_in_frame =
            contact_geometry_->findNearestPoint(transform_.shiftBaseStationToFrame(pnt), inside, normal);

        return transform_.shiftFrameStationToBase(closest_point_in_frame);
    }
    //=================================================================================================//
    GeometricShapeBrick::
        GeometricShapeBrick(const Vec3d &halfsize, SimTK::Transform transform, const std::string &shape_name)
        : GeometricShape(shape_name, transform), brick_(halfsize)
    {
        contact_geometry_ = &brick_;
    }
    //=================================================================================================//
    bool GeometricShapeBrick::checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED)
    {
        return brick_.getGeoBox().containsPoint(transform_.shiftBaseStationToFrame(pnt));
        ;
    }
    //=================================================================================================//
    Vec3d GeometricShapeBrick::findClosestPoint(const Vec3d &pnt)
    {
        bool inside = false;
        Vec3d closest_point_in_frame =
            brick_.getGeoBox().findClosestPointOnSurface(transform_.shiftBaseStationToFrame(pnt), inside);
        return transform_.shiftFrameStationToBase(closest_point_in_frame);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeBrick::findBounds()
    {
        Vec3d halfsize = brick_.getHalfLengths();

        //initial reference values
        Vec3d lower_bound = Vec3d(Infinity);
        Vec3d upper_bound = Vec3d(-Infinity);

        //lower left corner
        Vec3d lower_left = transform_.shiftFrameStationToBase(-halfsize);
        for (int j = 0; j != 3; ++j)
        {
            lower_bound[j] = SMIN(lower_bound[j], lower_left[j]);
            upper_bound[j] = SMAX(upper_bound[j], lower_left[j]);
        }

        //upper right corner
        Vec3d upper_right = transform_.shiftFrameStationToBase(halfsize);
        for (int j = 0; j != 3; ++j)
        {
            lower_bound[j] = SMIN(lower_bound[j], upper_right[j]);
            upper_bound[j] = SMAX(upper_bound[j], upper_right[j]);
        }

        return BoundingBox(lower_bound, upper_bound);
    }
    //=================================================================================================//
}