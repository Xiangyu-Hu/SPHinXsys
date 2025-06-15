#include "kernel_tabulated_ck.h"

namespace SPH
{
//=================================================================================================//
KernelTabulatedCK::KernelTabulatedCK(Kernel &kernel)
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

    dq_ = kernel.KernelSize() / Real(kernel_resolution_);
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

    delta_q_0_ = (-1.0 * dq_) * (-2.0 * dq_) * (-3.0 * dq_);
    delta_q_1_ = dq_ * (-1.0 * dq_) * (-2.0 * dq_);
    delta_q_2_ = (2.0 * dq_) * dq_ * (-1.0 * dq_);
    delta_q_3_ = (3.0 * dq_) * (2.0 * dq_) * dq_;
}
//=================================================================================================//
} // namespace SPH
