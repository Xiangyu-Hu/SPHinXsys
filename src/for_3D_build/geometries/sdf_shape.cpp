#include "sdf_shape.h"

namespace SPH
{
//=================================================================================================//
SDFShape::SDFShape(Real finest_grid_spacing, const std::string &shape_name)
    : Shape(shape_name), finest_grid_spacing_(finest_grid_spacing) {}
//=================================================================================================//
bool SDFShape::checkContain(const Vec3d &probe_point, bool BOUNDARY_INCLUDED)
{
    return probeSignedDistance(probe_point) <= (BOUNDARY_INCLUDED ? 0.0 : -SqrtEps);
}
//=================================================================================================//
Vec3d SDFShape::findClosestPoint(const Vec3d &probe_point)
{
    Real signed_distance = probeSignedDistance(probe_point);
    Vecd normal_direction = probeNormalDirection(probe_point);
    return probe_point - signed_distance * normal_direction;
}
//=================================================================================================//
Real SDFShape::probeSignedDistance(const Vec3d &probe_point)
{
    Real signed_distance = MaxReal;
    for (const auto &primitive_and_op : primitives_and_ops_)
    {
        SDFBase *sdf_entity = primitive_and_op.first;
        GeometricOps op = primitive_and_op.second;
        Real primitive_distance = (*sdf_entity)(probe_point);
        if (op == GeometricOps::add)
            signed_distance = SMIN(signed_distance, primitive_distance);
        else if (op == GeometricOps::sub)
            signed_distance = SMAX(signed_distance, -primitive_distance);
        else if (op == GeometricOps::sym_diff)
            signed_distance = ABS(primitive_distance);
        else if (op == GeometricOps::intersect)
            signed_distance = SMAX(signed_distance, primitive_distance);
    }
    return signed_distance;
}
//=================================================================================================//
Vecd SDFShape::probeNormalDirection(const Vec3d &probe_point)
{
    Real eps = 0.01 * finest_grid_spacing_;
    Vecd normal_direction;
    for (int i = 0; i != Dimensions; ++i)
    {
        Vecd probe_plus = probe_point;
        Vecd probe_minus = probe_point;
        probe_plus[i] += eps;
        probe_minus[i] -= eps;
        normal_direction[i] = (probeSignedDistance(probe_plus) -
                               probeSignedDistance(probe_minus)) /
                              (2.0 * eps);
    }
    return normal_direction.normalized();
}
//=================================================================================================//
BoundingBoxd SDFShape::findBounds() // only add and intersect operations are considered.
{
    BoundingBoxd bounding_box;
    for (const auto &primitive_and_op : primitives_and_ops_)
    {
        SDFBase *sdf_entity = primitive_and_op.first;
        BoundingBoxd primitive_bounds = sdf_entity->findBounds();
        GeometricOps op = primitive_and_op.second;
        if (op == GeometricOps::add)
            bounding_box = bounding_box.add(primitive_bounds);
        else if (op == GeometricOps::intersect)
            bounding_box = bounding_box.intersect(primitive_bounds);
    }

    return bounding_box;
}
//=================================================================================================//
} // namespace SPH
