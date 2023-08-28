#include "kernel_wenland_c2.h"

#include <cmath>

namespace SPH
{
//=================================================================================================//
KernelWendlandC2::KernelWendlandC2(Real h)
    : Kernel(h, 2.0, 2.0, "Wendland2CKernel")
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
DeviceKernelWendlandC2::DeviceKernelWendlandC2(Kernel& kernel)
{
    inv_h_ = 1.0 / kernel.SmoothingLength();
    factor_W_1D_ = kernel.FactorW1D();
    factor_W_2D_ = kernel.FactorW2D();
    factor_W_3D_ = kernel.FactorW3D();
    factor_dW_1D_ = inv_h_ * factor_W_1D_;
    factor_dW_2D_ = inv_h_ * factor_W_2D_;
    factor_dW_3D_ = inv_h_ * factor_W_3D_;
    rc_ref_ = kernel.CutOffRadius();
    rc_ref_sqr_ = kernel.CutOffRadiusSqr();
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::W(const DeviceReal &r_ij, const DeviceReal &displacement) const
{
    DeviceReal q = r_ij * inv_h_;
    return factor_W_1D_ * W_1D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::W(const DeviceReal &r_ij, const DeviceVec2d &displacement) const
{
    DeviceReal q = r_ij * inv_h_;
    return factor_W_2D_ * W_2D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::W(const DeviceReal &r_ij, const DeviceVec3d &displacement) const
{
    DeviceReal q = r_ij * inv_h_;
    return factor_W_3D_ * W_3D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::W_1D(const DeviceReal q) const
{
    return sycl::pow(1.0f - 0.5f * q, 4.0f) * (1.0f + 2.0f * q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::W_2D(const DeviceReal q) const
{
    return W_1D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::W_3D(const DeviceReal q) const
{
    return W_2D(q);
}

//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::dW(const DeviceReal &r_ij, const DeviceReal &displacement) const
{
    DeviceReal q = r_ij * inv_h_;
    return factor_dW_1D_ * dW_1D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::dW(const DeviceReal &r_ij, const DeviceVec2d &displacement) const
{
    DeviceReal q = r_ij * inv_h_;
    return factor_dW_2D_ * dW_2D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::dW(const DeviceReal &r_ij, const DeviceVec3d &displacement) const
{
    DeviceReal q = r_ij * inv_h_;
    return factor_dW_3D_ * dW_3D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::dW_1D(const DeviceReal q) const
{
    return 0.625f * sycl::pow(q - 2.0f, 3.0f) * q;
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::dW_2D(const DeviceReal q) const
{
    return dW_1D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::dW_3D(const DeviceReal q) const
{
    return dW_2D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::d2W_1D(const DeviceReal q) const
{
    return 1.25f * sycl::pow(q - 2.0f, 2.0f) * (2.0f * q - 1.0f);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::d2W_2D(const DeviceReal q) const
{
    return d2W_1D(q);
}
//=================================================================================================//
DeviceReal DeviceKernelWendlandC2::d2W_3D(const DeviceReal q) const
{
    return d2W_2D(q);
}
//=================================================================================================//
} // namespace SPH
