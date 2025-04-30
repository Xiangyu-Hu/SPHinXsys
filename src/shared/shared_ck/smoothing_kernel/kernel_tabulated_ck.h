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
 * @file kernel_tabulated_ck.h
 * @brief This is the class for tabulated kernel
 * which is applicable for all SPH kernels.
 * @author	Xiangyu Hu
 */

#ifndef KERNEL_TABULATED_CK_H
#define KERNEL_TABULATED_CK_H

#include "base_kernel.h"

namespace SPH
{
constexpr int KernelResolution = 20;
constexpr int TabulatedArraySize = KernelResolution + 4;

class KernelTabulatedCK
{
  public:
    explicit KernelTabulatedCK(Kernel &kernel);

    Real interpolateCubic(const Real *data, Real q) const
    {
        int location = (int)floor(q / dq_);
        int i = location + 1;
        Real shift1 = q - Real(location) * dq_; // shift1 correspond to i
        Real shift0 = shift1 + dq_;             // shift0 correspond to i-1
        Real shift2 = shift1 - dq_;             // shift2 correspond to i+1
        Real shift3 = shift1 - 2 * dq_;         ////shift3 correspond to i+2

        return (shift1 * shift2 * shift3) / dq0_ * data[i - 1] +
               (shift0 * shift2 * shift3) / dq1_ * data[i] +
               (shift0 * shift1 * shift3) / dq2_ * data[i + 1] +
               (shift0 * shift1 * shift2) / dq3_ * data[i + 2];
    };

    Real KernelSize() const { return kernel_size_; };

    Real W(const Real &displacement) const
    {
        Real q = displacement;
        return factor1d_ * interpolateCubic(w_1d, q);
    };

    Real W(const Vec2d &displacement) const
    {
        Real q = displacement.norm();
        return factor2d_ * interpolateCubic(w_1d, q);
    };

    Real W(const Vec3d &displacement) const
    {
        Real q = displacement.norm();
        return factor3d_ * interpolateCubic(w_1d, q);
    };

    Real dW(const Real &displacement) const
    {
        Real q = displacement;
        return factor1d_ * interpolateCubic(dw_1d, q);
    };
    Real dW(const Vec2d &displacement) const
    {
        Real q = displacement.norm();
        return factor2d_ * interpolateCubic(dw_1d, q);
    };
    Real dW(const Vec3d &displacement) const
    {
        Real q = displacement.norm();
        return factor3d_ * interpolateCubic(dw_1d, q);
    };

  private:
    Real kernel_size_;
    Real factor1d_, factor2d_, factor3d_;
    Real dq_, dq0_, dq1_, dq2_, dq3_;
    Real w_1d[TabulatedArraySize], dw_1d[TabulatedArraySize];
};
} // namespace SPH
#endif // KERNEL_TABULATED_CK_H