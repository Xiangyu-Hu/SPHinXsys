#include "mapping_shape.h"

namespace SPH
{
//=================================================================================================//
ExtrudeShape::ExtrudeShape(Shape &base_shape, Real thickness, const std::string &shape_name)
    : Shape(shape_name),
      thickness_(thickness), thickness_sqr_(thickness * thickness),
      base_shape_(base_shape) {};
//=================================================================================================//
Vecd ExtrudeShape::getShift(const Vecd &probe_point, const Vecd &original_closest_point)
{
    Vecd displacement = original_closest_point - probe_point;
    return thickness_ * displacement / (displacement.norm() + Eps);
}
//=================================================================================================//
bool ExtrudeShape::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    Vecd original_closest_point = base_shape_.findClosestPoint(probe_point);
    Vecd displacement = original_closest_point - probe_point;
    if (base_shape_.checkContain(probe_point))
    {
        return thickness_ > 0.0 ? true : displacement.squaredNorm() > thickness_sqr_;
    }
    else
    {
        return thickness_ < 0.0 ? false : displacement.squaredNorm() < thickness_sqr_;
    }
}
//=================================================================================================//
Vecd ExtrudeShape::findClosestPoint(const Vecd &probe_point)
{
    Vecd closest_point = base_shape_.findClosestPoint(probe_point);
    Vecd shift = getShift(probe_point, closest_point);
    closest_point += base_shape_.checkContain(probe_point) ? shift : -shift;
    return closest_point;
}
//=================================================================================================//
BoundingBoxd ExtrudeShape::findBounds()
{
    BoundingBoxd bounds = base_shape_.findBounds();
    return bounds.expand(thickness_ * Vecd::Ones());
}
//=================================================================================================//
} // namespace SPH
