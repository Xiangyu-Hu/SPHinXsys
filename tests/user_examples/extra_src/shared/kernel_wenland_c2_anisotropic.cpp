/**
 * @file 	kernel_wenland.cpp
 * @author	Luhui Han, Chi Zhang, Yongchuan Yu and Xiangyu Hu
 */
#include "kernel_wenland_c2_anisotropic.h"

#include <cmath>

namespace SPH
{

namespace Anisotropic
{
//=================================================================================================//
KernelWendlandC2::KernelWendlandC2(Real h)
    : Anisotropic::Kernel(h, 2.0, 2.0, "Wendland2CKernel")
{
    factor_W_1D_ = inv_h_ * 3.0 / 4.0;
    factor_W_2D_ = inv_h_ * inv_h_ * 7.0 / (4.0 * Pi);
    factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ * 21.0 / (16.0 * Pi);
    setDerivativeParameters();
}
//=================================================================================================//
Real KernelWendlandC2::W_1D(const Real q) const
{
    return pow(1.0 - 0.5 * q, 4) * (1.0 + 2.0 * q);
}
//=================================================================================================//
Real KernelWendlandC2::W_2D(const Real q) const
{
    return W_1D(q);
}
//=================================================================================================//
Real KernelWendlandC2::W_3D(const Real q) const
{
    return W_2D(q);
}
//=================================================================================================//
Real KernelWendlandC2::dW_1D(const Real q) const
{
    return 0.625 * pow(q - 2.0, 3) * q;
}
//=================================================================================================//
Real KernelWendlandC2::dW_2D(const Real q) const
{
    return dW_1D(q);
}
//=================================================================================================//
Real KernelWendlandC2::dW_3D(const Real q) const
{
    return dW_2D(q);
}
//=================================================================================================//
Real KernelWendlandC2::d2W_1D(const Real q) const
{
    return 1.25 * pow(q - 2.0, 2) * (2.0 * q - 1.0);
}
//=================================================================================================//
Real KernelWendlandC2::d2W_2D(const Real q) const
{
    return d2W_1D(q);
}
//=================================================================================================//
Real KernelWendlandC2::d2W_3D(const Real q) const
{
    return d2W_2D(q);
}
//=================================================================================================//
} // namespace Anisotropic
} // namespace SPH
