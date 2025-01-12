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
bool GeometricBox::checkContain(const Vecd &probe_point)
{
    bool is_contained = true;
    for (int i = 0; i != Dimensions; ++i)
    {
        if (ABS(probe_point[i]) > halfsize_[i]) // outside the box
        {
            is_contained = false;
            break;
        }
    }
    return is_contained;
}
//=================================================================================================//
Vecd GeometricBox::findClosestPoint(const Vecd &probe_point)
{
    Vecd c(probe_point);
    for (int i = 0; i < Dimensions; ++i)
    {
        if (c[i] < -halfsize_[i])
        {
            c[i] = -halfsize_[i];
        }
        else if (c[i] > halfsize_[i])
        {
            c[i] = halfsize_[i];
        }
    }
    return c;
}
//=================================================================================================//
} // namespace SPH