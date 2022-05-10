/**
 * @file 	geometric_shape.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "geometric_shape.h"

namespace SPH
{
    //=================================================================================================//
    GeometricShapeCircle::GeometricShapeCircle(const Vec2d &center, Real radius,
                                               const std::string &shape_name)
        : Shape(shape_name), center_(center), radius_(radius) {}
    //=================================================================================================//
    bool GeometricShapeCircle::checkContain(const Vec2d &input_pnt, bool BOUNDARY_INCLUDED)
    {
        return (input_pnt - center_).norm() < radius_;
    }
    //=================================================================================================//
    Vec2d GeometricShapeCircle::findClosestPoint(const Vec2d &input_pnt)
    {
        Vec2d displacement = input_pnt - center_;
        Real distance = displacement.norm();
        Real cosine = (SGN(displacement[0])*(ABS(displacement[0])) + TinyReal) / (distance + TinyReal);
        Real sine = displacement[1] / (distance + TinyReal);
        return input_pnt + (radius_ - distance) * Vec2d(cosine, sine);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeCircle::findBounds()
    {
        Vec2d shift = Vec2d(radius_, radius_);
        return BoundingBox(center_ - shift, center_ + shift);
    }
    //=================================================================================================//
    Real GeometricShapeCircle::findSignedDistance(const Vecd &input_pnt)
    {
        return (input_pnt - center_).norm() - radius_;
    }
    //=================================================================================================//
}