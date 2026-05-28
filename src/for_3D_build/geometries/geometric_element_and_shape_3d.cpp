#include "geometric_element_and_shape_3d.h"
#include "geometric_shape.h"

#include "io_environment.h"
#include "triangle_mesh_shape.h"

namespace SPH
{
//=================================================================================================//
GeometricCylinder::GeometricCylinder(Real radius, Real halflength)
    : radius_(radius), halflength_(halflength)
{
    if (radius < 0.0)
    {
        std::cout << "\n Error: the GeometricCylinder radius must be positive! " << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    if (halflength < 0.0)
    {
        std::cout << "\n Error: the GeometricCylinder halflength must be positive! " << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
Vecd GeometricCylinder::findClosestPoint(const Vecd &probe_point)
{
    Real axial = probe_point[0];
    Real radial_distance = probe_point.tail(Dimensions - 1).norm();

    Real dh = ABS(axial) - halflength_;
    Real dr = radial_distance - radius_;
    Real clamped_axial = clamp(axial, -halflength_, halflength_);

    Vecd result;
    result[0] = clamped_axial;
    if (radial_distance > Eps)
        result.tail(Dimensions - 1) = (radius_ / radial_distance) * probe_point.tail(Dimensions - 1);
    else
    {
        result.tail(Dimensions - 1).setZero();
        result.tail(Dimensions - 1)[0] = radius_;
    }

    if (SMAX(dr, dh) <= 0.0)
    {
        // Point inside: project to nearest surface
        if (dr >= dh)
            result[0] = axial; // closer to cylindrical surface
        else if (radial_distance > Eps)
            result.tail(Dimensions - 1) = probe_point.tail(Dimensions - 1); // closer to end cap
        else
            result.tail(Dimensions - 1).setZero();
    }
    else if (dh > 0.0 && dr <= 0.0)
    {
        // Outside axially, inside radially: project to cap
        if (radial_distance > Eps)
            result.tail(Dimensions - 1) = probe_point.tail(Dimensions - 1);
        else
            result.tail(Dimensions - 1).setZero();
    }
    else if (dr > 0.0 && dh <= 0.0)
    {
        // Outside radially, inside axially: project to cylindrical surface
        result[0] = axial;
    }
    // else: outside both, result already holds (clamped_axial, normalized_radial)

    return result;
}
//=================================================================================================//
BoundingBox3d GeometricCylinder::findBounds()
{
    Vecd min_corner = Vecd::Constant(-radius_);
    Vecd max_corner = Vecd::Constant(radius_);
    min_corner[0] = -halflength_;
    max_corner[0] = halflength_;
    return BoundingBox3d(min_corner, max_corner);
}
//=================================================================================================//
GeometricShapeCylinder::GeometricShapeCylinder(const Transform &transform, Real radius, Real halflength,
                                               const std::string &name)
    : TransformShape<GeometricCylinder>(name, transform, radius, halflength) {}
//=================================================================================================//
void GeometricShapeBox::writeGeometricShapeBoxToVtp(Real scale_factor)
{
    TriangleMeshShapeBrick shape(HalfSize(), 1, Vecd::Zero(), Name());
    shape.writTriangleMeshShapeToVtp(getTransform(), scale_factor);
}
//=================================================================================================//
} // namespace SPH