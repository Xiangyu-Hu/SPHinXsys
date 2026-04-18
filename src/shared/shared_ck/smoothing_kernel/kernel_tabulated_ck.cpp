#include "kernel_tabulated_ck.h"

namespace SPH
{
//=================================================================================================//
KernelTabulatedCK::KernelTabulatedCK(Kernel &kernel)
{
    kernel_size_ = kernel.KernelSize();
    dimension_factor_1D_ = kernel.DimensionFactor1D();
    dimension_factor_2D_ = kernel.DimensionFactor2D();
    dimension_factor_3D_ = kernel.DimensionFactor3D();

    dq_ = kernel_size_ / Real(kernel_resolution_);
    for (int i = 0; i < tabulated_size_; i++)
    {
        w_1d[i] = kernel.W_1D(Real(i - 1) * dq_);
        dw_1d[i] = kernel.dW_1D(Real(i - 1) * dq_);
        d2w_1d[i] = kernel.d2W_1D(Real(i - 1) * dq_);
    }

    delta_q_0_ = (-1.0 * dq_) * (-2.0 * dq_) * (-3.0 * dq_);
    delta_q_1_ = dq_ * (-1.0 * dq_) * (-2.0 * dq_);
    delta_q_2_ = (2.0 * dq_) * dq_ * (-1.0 * dq_);
    delta_q_3_ = (3.0 * dq_) * (2.0 * dq_) * dq_;
}
//=================================================================================================//
} // namespace SPH
