#include "kernel_tabulated_ck.h"

namespace SPH
{
//=================================================================================================//
KernelTabulatedCK::KernelTabulatedCK(Kernel &kernel)
    : kernel_size_(kernel.KernelSize()),
      factor1d_(kernel.FactorW1D()), factor2d_(kernel.FactorW2D()), factor3d_(kernel.FactorW3D()),
      dq_(kernel_size_ / Real(KernelResolution)),
      dq0_((-1.0 * dq_) * (-2.0 * dq_) * (-3.0 * dq_)),
      dq1_(dq_ * (-1.0 * dq_) * (-2.0 * dq_)),
      dq2_((2.0 * dq_) * dq_ * (-1.0 * dq_)),
      dq3_((3.0 * dq_) * (2.0 * dq_) * dq_)
{
    for (int i = 0; i < KernelResolution; i++)
    {
        w_1d[i] = kernel.W_1D(Real(i - 1) * dq_);
        dw_1d[i] = kernel.dW_1D(Real(i - 1) * dq_);
    }
    // kernel trailing zeros
    for (int i = KernelResolution; i < TabulatedArraySize; i++)
    {
        w_1d[i] = 0.0;
        dw_1d[i] = 0.0;
    }
} //=================================================================================================//
} // namespace SPH
