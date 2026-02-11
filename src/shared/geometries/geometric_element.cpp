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
    Real axial_projection = probe_point.dot(axis_);
    Vecd radial_vector = probe_point - axial_projection * axis_;
    Real radial_distance = radial_vector.norm();
    
    // Determine if the point is within the axial extent
    bool within_axial_extent = ABS(axial_projection) <= halflength_;
    // Determine if the point is within the radial extent
    bool within_radial_extent = radial_distance <= radius_;
    
    if (within_axial_extent && within_radial_extent)
    {
        // Point is inside the cylinder, find closest surface point
        Real distance_to_side = radius_ - radial_distance;
        Real distance_to_cap = halflength_ - ABS(axial_projection);
        
        if (distance_to_side < distance_to_cap)
        {
            // Closest to the cylindrical surface
            if (radial_distance > Eps)
            {
                return axial_projection * axis_ + radius_ * radial_vector / radial_distance;
            }
            else
            {
                // Point is on the axis, choose any radial direction
                Vecd radial_out = Vecd::Zero();
                radial_out[0] = radius_;
                return axial_projection * axis_ + radial_out;
            }
        }
        else
        {
            // Closest to one of the end caps
            Real cap_position = axial_projection > 0 ? halflength_ : -halflength_;
            return cap_position * axis_ + radial_vector;
        }
    }
    else if (within_axial_extent)
    {
        // Outside radially, but within axial extent
        if (radial_distance > Eps)
        {
            return axial_projection * axis_ + radius_ * radial_vector / radial_distance;
        }
        else
        {
            // Point is on the axis but outside
            Vecd radial_out = Vecd::Zero();
            radial_out[0] = radius_;
            return axial_projection * axis_ + radial_out;
        }
    }
    else if (within_radial_extent)
    {
        // Outside axially, but within radial extent
        Real cap_position = axial_projection > 0 ? halflength_ : -halflength_;
        return cap_position * axis_ + radial_vector;
    }
    else
    {
        // Outside both axially and radially - closest to edge of end cap
        Real cap_position = axial_projection > 0 ? halflength_ : -halflength_;
        if (radial_distance > Eps)
        {
            return cap_position * axis_ + radius_ * radial_vector / radial_distance;
        }
        else
        {
            Vecd radial_out = Vecd::Zero();
            radial_out[0] = radius_;
            return cap_position * axis_ + radial_out;
        }
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