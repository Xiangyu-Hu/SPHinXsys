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
    return radial_vector.norm() <= radius_;
}
//=================================================================================================//
Vecd GeometricCylinder::findClosestPoint(const Vecd &probe_point)
{
    // Decompose probe point into axial and radial components
    Real axial_projection = probe_point.dot(axis_);
    Vecd radial_vector = probe_point - axial_projection * axis_;
    Real radial_distance = radial_vector.norm();
    
    // Compute signed distance components
    // dh: signed distance along the axis (positive if outside caps)
    Real dh = ABS(axial_projection) - halflength_;
    // dr: signed distance normal to the axis (positive if outside cylinder surface)
    Real dr = radial_distance - radius_;
    
    // Clamp axial projection to cylinder caps using clamp function
    Real clamped_axial = clamp(axial_projection, -halflength_, halflength_);
    
    // Normalize radial vector and scale to radius
    Vecd normalized_radial;
    if (radial_distance > Eps)
        normalized_radial = radius_ * radial_vector / radial_distance;
    else
    {
        // Point on axis - choose arbitrary radial direction
        normalized_radial = Vecd::Zero();
        normalized_radial[0] = radius_;
    }
    
    // Determine closest point based on signed distances using SMAX
    if (SMAX(dr, dh) <= 0.0)
    {
        // Point inside cylinder - project to nearest surface
        if (SMAX(dr, dh) == dr)
        {
            // Closer to cylindrical surface
            return axial_projection * axis_ + normalized_radial;
        }
        else
        {
            // Closer to end cap
            return clamped_axial * axis_ + radial_vector;
        }
    }
    else if (dr > 0.0 && dh <= 0.0)
    {
        // Outside radially, inside axially
        return axial_projection * axis_ + normalized_radial;
    }
    else if (dr <= 0.0 && dh > 0.0)
    {
        // Inside radially, outside axially
        return clamped_axial * axis_ + radial_vector;
    }
    else
    {
        // Outside both axially and radially
        return clamped_axial * axis_ + normalized_radial;
    }
}
//=================================================================================================//
BoundingBoxd GeometricCylinder::findBounds()
{
    // Create bounding box that encompasses the cylinder
    // For each dimension i, the extent is:
    // halflength * |axis[i]| + radius * sqrt(1 - axis[i]^2)
    // This accounts for the cylinder extending halflength along the axis
    // and radius in directions perpendicular to the axis
    
    Vecd min_corner = Vecd::Zero();
    Vecd max_corner = Vecd::Zero();
    
    for (int i = 0; i < Dimensions; ++i)
    {
        Real axis_component = ABS(axis_[i]);
        Real perpendicular_factor = std::sqrt(1.0 - axis_component * axis_component);
        Real extent = halflength_ * axis_component + radius_ * perpendicular_factor;
        min_corner[i] = -extent;
        max_corner[i] = extent;
    }
    
    return BoundingBoxd(min_corner, max_corner);
}
//=================================================================================================//
} // namespace SPH