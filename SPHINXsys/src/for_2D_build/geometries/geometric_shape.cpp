/**
 * @file 	geometric_shape.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "geometric_shape.h"

namespace SPH
{
    //=================================================================================================//
    GeometricShapeBox::GeometricShapeBox(const Vec2d &halfsize, const std::string &shape_name)
        : Shape(shape_name), halfsize_(halfsize),
          multi_polygon_({-halfsize, Vec2d(-halfsize[0], halfsize[1]), halfsize,
                          Vec2d(halfsize[0], -halfsize[1]), -halfsize})
    {
        if (halfsize[0] < 0.0 || halfsize[1] < 0.0)
        {
            std::cout << "\n Error: the GeometricShapeBox half size must be positive! " << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }
    //=================================================================================================//
    bool GeometricShapeBox::checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED)
    {
        return ABS(pnt[0]) < halfsize_[0] && ABS(pnt[1]) < halfsize_[1];
    }
    //=================================================================================================//
    Vec2d GeometricShapeBox::findClosestPoint(const Vec2d &pnt)
    {
        return multi_polygon_.findClosestPoint(pnt);
        ;
    }
    //=================================================================================================//
    BoundingBox GeometricShapeBox::findBounds()
    {
        return BoundingBox(-halfsize_, halfsize_);
    }
    //=================================================================================================//
    GeometricShapeBall::GeometricShapeBall(const Vec2d &center, Real radius,
                                           const std::string &shape_name)
        : Shape(shape_name), center_(center), radius_(radius) {}
    //=================================================================================================//
    bool GeometricShapeBall::checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED)
    {
        return (pnt - center_).norm() < radius_;
    }
    //=================================================================================================//
    Vec2d GeometricShapeBall::findClosestPoint(const Vec2d &pnt)
    {
        Vec2d displacement = pnt - center_;
        Real distance = displacement.norm();
        Real cosine = (SGN(displacement[0]) * (ABS(displacement[0])) + TinyReal) / (distance + TinyReal);
        Real sine = displacement[1] / (distance + TinyReal);
        return pnt + (radius_ - distance) * Vec2d(cosine, sine);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeBall::findBounds()
    {
        Vec2d shift = Vec2d(radius_, radius_);
        return BoundingBox(center_ - shift, center_ + shift);
    }
    //=================================================================================================//
}