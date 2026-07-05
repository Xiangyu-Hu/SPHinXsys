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
} // namespace SPH