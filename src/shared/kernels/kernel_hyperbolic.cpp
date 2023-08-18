#include "kernel_hyperbolic.h"

#include <cmath>

namespace SPH
{
//=================================================================================================//
KernelHyperbolic::KernelHyperbolic(Real h)
    : Kernel(h, 2.0, 2.0, "HyperbolicKernel")
{
    factor_W_1D_ = inv_h_ / 7.0;
    factor_W_2D_ = inv_h_ * inv_h_ / (3.0 * Pi);
    factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ * 15.0 / (62.0 * Pi);
    setDerivativeParameters();
}
//=================================================================================================//
Real KernelHyperbolic::W_1D(const Real q) const
{
    if (q < 1.0)
    {
        return (6.0 - 6.0 * q + pow(q, 3));
    }
    else
    {
        return pow(2.0 - q, 3);
    }
}
//=================================================================================================//
Real KernelHyperbolic::W_2D(const Real q) const
{
    return W_1D(q);
}
//=================================================================================================//
Real KernelHyperbolic::W_3D(const Real q) const
{
    return W_1D(q);
}
//=================================================================================================//
Real KernelHyperbolic::dW_1D(const Real q) const
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
Real KernelHyperbolic::dW_2D(const Real q) const
{
    return dW_1D(q);
}
//=================================================================================================//
Real KernelHyperbolic::dW_3D(const Real q) const
{
    return dW_1D(q);
}
//=================================================================================================//
Real KernelHyperbolic::d2W_1D(const Real q) const
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
Real KernelHyperbolic::d2W_2D(const Real q) const
{
    return d2W_1D(q);
}
//=================================================================================================//
Real KernelHyperbolic::d2W_3D(const Real q) const
{
    return d2W_1D(q);
}
//=================================================================================================//
} // namespace SPH
