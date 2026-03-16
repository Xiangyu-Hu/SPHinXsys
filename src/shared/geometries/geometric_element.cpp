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
BoundingBoxd GeometricCylinder::findBounds()
{
    Vecd min_corner = Vecd::Constant(-radius_);
    Vecd max_corner = Vecd::Constant(radius_);
    min_corner[0] = -halflength_;
    max_corner[0] = halflength_;
    return BoundingBoxd(min_corner, max_corner);
}
//=================================================================================================//
} // namespace SPH