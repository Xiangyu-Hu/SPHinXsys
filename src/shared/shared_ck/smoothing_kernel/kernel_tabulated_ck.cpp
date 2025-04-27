#include "kernel_tabulated_ck.h"

namespace SPH
{
//=================================================================================================//
KernelTabulatedCK::KernelTabulatedCK(Kernel &kernel)
    : kernel_size_(kernel.KernelSize()),
      factor1D_(kernel.FactorW1D()), factor2D_(kernel.FactorW2D()), factor3D_(kernel.FactorW3D()),
      dq_(kernel_size_ / Real(kernel_resolution_)),
      delta_q_0_((-1.0 * dq_) * (-2.0 * dq_) * (-3.0 * dq_)),
      delta_q_1_(dq_ * (-1.0 * dq_) * (-2.0 * dq_)),
      delta_q_2_((2.0 * dq_) * dq_ * (-1.0 * dq_)),
      delta_q_3_((3.0 * dq_) * (2.0 * dq_) * dq_)
{
    for (int i = 0; i < kernel_resolution_; i++)
    {
        w_1d[i] = kernel.W_1D(Real(i - 1) * dq_);
        dw_1d[i] = kernel.dW_1D(Real(i - 1) * dq_);
    }
    // kernel trailing zeros
    for (int i = kernel_resolution_; i < tabulated_array_size_; i++)
    {
        w_1d[i] = 0.0;
        dw_1d[i] = 0.0;
    }
} //=================================================================================================//
} // namespace SPH
