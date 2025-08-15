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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
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
constexpr int kernel_resolution_ = 20;
constexpr int tabulated_size_ = kernel_resolution_ + 4;

class KernelTabulatedCK
{
  public:
    explicit KernelTabulatedCK(Kernel &kernel);

    Real interpolateCubic(const Real *data, Real q) const
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

    Real factorW1D(const Real &h_ratio) const { return h_ratio; };
    Real factorW2D(const Real &h_ratio) const { return h_ratio * h_ratio; };
    Real factorW3D(const Real &h_ratio) const { return h_ratio * h_ratio * h_ratio; };
    Real factordW1D(const Real &h_ratio) const { return factorW1D(h_ratio) * h_ratio; };
    Real factordW2D(const Real &h_ratio) const { return factorW2D(h_ratio) * h_ratio; };
    Real factordW3D(const Real &h_ratio) const { return factorW3D(h_ratio) * h_ratio; };
    Real factord2W1D(const Real &h_ratio) const { return factordW1D(h_ratio) * h_ratio; };
    Real factord2W2D(const Real &h_ratio) const { return factordW2D(h_ratio) * h_ratio; };
    Real factord2W3D(const Real &h_ratio) const { return factordW3D(h_ratio) * h_ratio; };

    Real W(const Real &displacement) const
    {
        Real q = displacement * inv_h_;
        return factor_W_1D_ * interpolateCubic(w_1d, q);
    };

    Real W(const Vec2d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_W_2D_ * interpolateCubic(w_1d, q);
    };

    Real W(const Vec3d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_W_3D_ * interpolateCubic(w_1d, q);
    };

    Real W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_W_1D_ * interpolateCubic(w_1d, q) * factorW1D(h_ratio);
    }

    Real W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_W_2D_ * interpolateCubic(w_1d, q) * factorW2D(h_ratio);
    }

    Real W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_W_3D_ * interpolateCubic(w_1d, q) * factorW3D(h_ratio);
    }

    Real dW(const Real &displacement) const
    {
        Real q = displacement * inv_h_;
        return factor_dW_1D_ * interpolateCubic(dw_1d, q);
    };
    Real dW(const Vec2d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_dW_2D_ * interpolateCubic(dw_1d, q);
    };
    Real dW(const Vec3d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_dW_3D_ * interpolateCubic(dw_1d, q);
    };

    Real dW(const Real &h_ratio, const Real &r_ij, const Real &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_dW_1D_ * interpolateCubic(dw_1d, q) * factordW1D(h_ratio);
    };
    Real dW(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_dW_2D_ * interpolateCubic(dw_1d, q) * factordW2D(h_ratio);
    };
    Real dW(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_dW_3D_ * interpolateCubic(dw_1d, q) * factordW3D(h_ratio);
    };

    Real d2W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_d2W_1D_ * interpolateCubic(d2w_1d, q) * factord2W1D(h_ratio);
    }

    Real d2W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_d2W_2D_ * interpolateCubic(d2w_1d, q) * factord2W2D(h_ratio);
    }

    Real d2W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const
    {
        Real q = r_ij * inv_h_ * h_ratio;
        return factor_d2W_3D_ * interpolateCubic(d2w_1d, q) * factord2W3D(h_ratio);
    }

    Real CutOffRadius() const { return rc_ref_; };
    Real CutOffRadiusSqr() const { return rc_ref_sqr_; };
    Real KernelSize() const { return kernel_size_; };

    Real CutOffRadius(Real h_ratio) const { return rc_ref_ / h_ratio; };

    inline Real DimensionFactor(const Real &) const { return dimension_factor_1D_; };
    inline Real DimensionFactor(const Vec2d &) const { return dimension_factor_2D_; };
    inline Real DimensionFactor(const Vec3d &) const { return dimension_factor_3D_; };

    inline Real normalized_W(Real normalized_distance) const
    {
        return interpolateCubic(w_1d, normalized_distance);
    };

    inline Real normalized_dW(Real normalized_distance) const
    {
        return interpolateCubic(dw_1d, normalized_distance);
    };

    inline Real normalized_d2W(Real normalized_distance) const
    {
        return interpolateCubic(d2w_1d, normalized_distance);
    };

  protected:
    Real dimension_factor_1D_, dimension_factor_2D_, dimension_factor_3D_;

  private:
    Real inv_h_, rc_ref_, rc_ref_sqr_, kernel_size_;
    Real factor_W_1D_, factor_W_2D_, factor_W_3D_;
    Real factor_dW_1D_, factor_dW_2D_, factor_dW_3D_;
    Real factor_d2W_1D_, factor_d2W_2D_, factor_d2W_3D_;
    Real dq_, delta_q_0_, delta_q_1_, delta_q_2_, delta_q_3_;
    Real w_1d[tabulated_size_], dw_1d[tabulated_size_], d2w_1d[tabulated_size_];
};
} // namespace SPH
#endif // KERNEL_TABULATED_CK_H