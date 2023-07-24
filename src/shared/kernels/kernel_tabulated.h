/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	kernel_tabulated.h
 * @brief 	This is the class for tabulated kernels using template.
 * @author	Yongchuan Yu, Massoud Rezevand, Chi Zhang and Xiangyu Hu
 */

#ifndef KERNEL_TABULATED_H
#define KERNEL_TABULATED_H

#include "base_kernel.h"
#include <cmath>

namespace SPH
{
template <class KernelType>
class KernelTabulated : public Kernel
{
  protected:
    KernelType original_kernel_;
    int kernel_resolution_;
    Real dq_, delta_q_0_, delta_q_1_, delta_q_2_, delta_q_3_;
    StdVec<Real> w_1d, w_2d, w_3d;
    StdVec<Real> dw_1d, dw_2d, dw_3d;
    StdVec<Real> d2w_1d, d2w_2d, d2w_3d;

    /** interpolation function, Four-point Lagrangian interpolation. */
    Real InterpolationCubic(const StdVec<Real> &data, Real q) const
    {
        int location = (int)floor(q / dq_);
        int i = location + 1;
        Real fraction_1 = q - Real(location) * dq_; // fraction_1 correspond to i
        Real fraction_0 = fraction_1 + dq_;         // fraction_0 correspond to i-1
        Real fraction_2 = fraction_1 - dq_;         // fraction_2 correspond to i+1
        Real fraction_3 = fraction_1 - 2 * dq_;     ////fraction_3 correspond to i+2

        return (fraction_1 * fraction_2 * fraction_3) / delta_q_0_ * data[i - 1] +
               (fraction_0 * fraction_2 * fraction_3) / delta_q_1_ * data[i] +
               (fraction_0 * fraction_1 * fraction_3) / delta_q_2_ * data[i + 1] +
               (fraction_0 * fraction_1 * fraction_2) / delta_q_3_ * data[i + 2];
    };

  public:
    explicit KernelTabulated(Real h, int kernel_resolution);

    virtual Real W_1D(const Real q) const override;
    virtual Real W_2D(const Real q) const override;
    virtual Real W_3D(const Real q) const override;

    virtual Real dW_1D(const Real q) const override;
    virtual Real dW_2D(const Real q) const override;
    virtual Real dW_3D(const Real q) const override;

    virtual Real d2W_1D(const Real q) const override;
    virtual Real d2W_2D(const Real q) const override;
    virtual Real d2W_3D(const Real q) const override;
};
//=================================================================================================//
template <class KernelType>
KernelTabulated<KernelType>::KernelTabulated(Real h, int kernel_resolution)
    : Kernel(h, 2.0, 2.0, "Tabulated"), original_kernel_(h),
      kernel_resolution_(kernel_resolution)
{
    kernel_name_ += original_kernel_.Name();
    kernel_size_ = original_kernel_.KernelSize();
    truncation_ = original_kernel_.Truncation();
    rc_ref_ = original_kernel_.CutOffRadius();
    rc_ref_sqr_ = original_kernel_.CutOffRadiusSqr();

    factor_W_1D_ = original_kernel_.FactorW1D();
    factor_W_2D_ = original_kernel_.FactorW2D();
    factor_W_3D_ = original_kernel_.FactorW3D();

    setDerivativeParameters();

    dq_ = KernelSize() / Real(kernel_resolution_);
    for (int i = 0; i < kernel_resolution_ + 4; i++)
    {
        w_1d.push_back(original_kernel_.W_1D(Real(i - 1) * dq_));
        w_2d.push_back(original_kernel_.W_2D(Real(i - 1) * dq_));
        w_3d.push_back(original_kernel_.W_3D(Real(i - 1) * dq_));
        dw_1d.push_back(original_kernel_.dW_1D(Real(i - 1) * dq_));
        dw_2d.push_back(original_kernel_.dW_2D(Real(i - 1) * dq_));
        dw_3d.push_back(original_kernel_.dW_3D(Real(i - 1) * dq_));
        d2w_1d.push_back(original_kernel_.d2W_1D(Real(i - 1) * dq_));
        d2w_2d.push_back(original_kernel_.d2W_2D(Real(i - 1) * dq_));
        d2w_3d.push_back(original_kernel_.d2W_3D(Real(i - 1) * dq_));
    }

    delta_q_0_ = (-1.0 * dq_) * (-2.0 * dq_) * (-3.0 * dq_);
    delta_q_1_ = dq_ * (-1.0 * dq_) * (-2.0 * dq_);
    delta_q_2_ = (2.0 * dq_) * dq_ * (-1.0 * dq_);
    delta_q_3_ = (3.0 * dq_) * (2.0 * dq_) * dq_;
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::W_1D(Real q) const
{
    return InterpolationCubic(w_1d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::W_2D(Real q) const
{
    return InterpolationCubic(w_2d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::W_3D(Real q) const
{
    return InterpolationCubic(w_3d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::dW_1D(Real q) const
{
    return InterpolationCubic(dw_1d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::dW_2D(Real q) const
{
    return InterpolationCubic(dw_2d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::dW_3D(Real q) const
{
    return InterpolationCubic(dw_3d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::d2W_1D(Real q) const
{
    return InterpolationCubic(d2w_1d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::d2W_2D(Real q) const
{
    return InterpolationCubic(d2w_2d, q);
}
//=================================================================================================//
template <class KernelType>
Real KernelTabulated<KernelType>::d2W_3D(Real q) const
{
    return InterpolationCubic(d2w_3d, q);
}
//=================================================================================================//
} // namespace SPH
#endif // KERNEL_TABULATED_H