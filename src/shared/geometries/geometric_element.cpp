#include "geometric_element.h"

namespace SPH
{
//=================================================================================================//
GeometricBox::GeometricBox(const Vecd &halfsize) : halfsize_(halfsize)
{
    for (int i = 0; i != Dimensions; ++i)
    {
        if (halfsize[i] < 0.0)
        {
            std::cout << "\n Error: the GeometricBox half size must be positive! " << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }
}
//=================================================================================================//
Vecd GeometricBox::findClosestPoint(const Vecd &probe_point)
{
    // This is from Simbody, Geo_Box.h
    //  first step to find the closest point for the point outside the box
    Vecd c(probe_point); // tentatively inside
    bool ptWasInside = true;
    for (int i = 0; i < Dimensions; ++i)
    {
        if (c[i] < -halfsize_[i])
        {
            c[i] = -halfsize_[i];
            ptWasInside = false;
        }
        else if (c[i] > halfsize_[i])
        {
            c[i] = halfsize_[i];
            ptWasInside = false;
        }
    }

    // second step to find the closest point for the point inside the box
    if (ptWasInside)
    {
        Real dToSide = halfsize_[0] - std::abs(c[0]);

        int which = 0; // tentatively the nearest side
        Real minDist = dToSide;
        for (int i = 1; i < Dimensions; ++i)
        {
            dToSide = halfsize_[i] - std::abs(c[i]);
            if (dToSide < minDist)
            {
                which = i;
                minDist = dToSide;
            }
        }
        // Now project the point to the nearest side.
        c[which] = c[which] < 0 ? -halfsize_[which] : halfsize_[which];
    }
    return c;
}
//=================================================================================================//
BoundingBoxd GeometricBox::findBounds()
{
    return BoundingBoxd(-halfsize_, halfsize_);
}
//=================================================================================================//
GeometricBall::GeometricBall(Real radius) : radius_(radius)
{
    if (radius < 0.0)
    {
        std::cout << "\n Error: the GeometricBall radius must be positive! " << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
bool GeometricBall::checkContain(const Vecd &probe_point)
{
    return probe_point.norm() < radius_;
}
//=================================================================================================//
Vecd GeometricBall::findClosestPoint(const Vecd &probe_point)
{
    return radius_ * probe_point.normalized();
}
//=================================================================================================//
BoundingBoxd GeometricBall::findBounds()
{
    Vecd shift = radius_ * Vecd::Ones();
    return BoundingBoxd(-shift, shift);
}
//=================================================================================================//
GeometricCylinder::GeometricCylinder(const Vecd &axis, Real radius, Real halflength)
    : axis_(axis.normalized()), radius_(radius), halflength_(halflength)
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
bool GeometricCylinder::checkContain(const Vecd &probe_point)
{
    Real axial_projection = probe_point.dot(axis_);
    if (ABS(axial_projection) > halflength_)
    {
        return false;
    }
    Vecd radial_vector = probe_point - axial_projection * axis_;
    return radial_vector.norm() < radius_;
}
//=================================================================================================//
Vecd GeometricCylinder::findClosestPoint(const Vecd &probe_point)
{
    Real axial_projection = probe_point.dot(axis_);
    Vecd radial_vector = probe_point - axial_projection * axis_;
    Real radial_distance = radial_vector.norm();
    
    // Clamp axial projection to cylinder length
    Real clamped_axial = axial_projection;
    if (clamped_axial < -halflength_)
        clamped_axial = -halflength_;
    else if (clamped_axial > halflength_)
        clamped_axial = halflength_;
    
    // Project radial component to cylinder radius
    Vecd clamped_radial = radial_vector;
    if (radial_distance > Eps)
    {
        clamped_radial = radius_ * radial_vector / radial_distance;
    }
    else
    {
        // Point is on the axis, choose any radial direction
        clamped_radial = Vecd::Zero();
        clamped_radial[0] = radius_;
    }
    
    return clamped_axial * axis_ + clamped_radial;
}
//=================================================================================================//
BoundingBoxd GeometricCylinder::findBounds()
{
    // Create bounding box that encompasses the cylinder
    // This is a conservative approximation
    Vecd half_axis = halflength_ * axis_;
    Vecd radial_extent = radius_ * Vecd::Ones();
    
    Vecd min_corner = -ABS(half_axis[0]) * Vecd::Unit(0);
    Vecd max_corner = ABS(half_axis[0]) * Vecd::Unit(0);
    
    for (int i = 0; i < Dimensions; ++i)
    {
        min_corner[i] = -ABS(half_axis[i]) - radius_;
        max_corner[i] = ABS(half_axis[i]) + radius_;
    }
    
    return BoundingBoxd(min_corner, max_corner);
}
//=================================================================================================//
} // namespace SPH