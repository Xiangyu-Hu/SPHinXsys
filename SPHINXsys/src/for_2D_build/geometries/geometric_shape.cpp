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
        Real level_set = distance - radius_ ;
        Real cosine = displacement[0] / (distance + TinyReal);
        Real sine_abs = sqrt(1.0 - cosine * cosine);
        Real sine = displacement[1] > 0.0 ?  sine_abs : -sine_abs;
       return input_pnt - level_set * Vec2d(cosine, sine);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeCircle::findBounds()
    {
        Vec2d shift = Vec2d(0.70710678118, 0.70710678118) * radius_;
        return BoundingBox(center_ - shift,  center_ + shift);
    }
    //=================================================================================================//
}