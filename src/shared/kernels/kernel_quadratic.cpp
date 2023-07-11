#include "kernel_quadratic.h"

#include <cmath>

namespace SPH
{
//=================================================================================================//
KernelQuadratic::KernelQuadratic(Real h)
    : Kernel(h, 2.0, 2.0, "QuadraticKernel")
{
    factor_W_1D_ = inv_h_ / 7.0;
    factor_W_2D_ = inv_h_ * inv_h_ / (3.0 * Pi);
    factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ / Pi;
    setDerivativeParameters();
}
//=================================================================================================//
Real KernelQuadratic::W_1D(const Real q) const
{
    return 5.0 * (3.0 * q * q - 12.0 * q + 12.0) / 64.0;
}
//=================================================================================================//
Real KernelQuadratic::W_2D(const Real q) const
{
    return W_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::W_3D(const Real q) const
{
    return W_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::dW_1D(const Real q) const
{
    if (q < 1.0)
    {
        return (-6.0 + 3.0 * pow(q, 2));
    }
    else
    {
        return pow(2.0 - q, 2) * (-1.0);
    }
}
//=================================================================================================//
Real KernelQuadratic::dW_2D(const Real q) const
{
    return dW_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::dW_3D(const Real q) const
{
    return 15.0 * (q - 2.0) / 32.0;
}
//=================================================================================================//
Real KernelQuadratic::d2W_1D(const Real q) const
{
    if (q < 1.0)
    {
        return 6.0 * q;
    }
    else
    {
        return 2.0 * (2.0 - q);
    }
}
//=================================================================================================//
Real KernelQuadratic::d2W_2D(const Real q) const
{
    return d2W_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::d2W_3D(const Real q) const
{
    return 15.0 / 32.0;
}
//=================================================================================================//
} // namespace SPH
