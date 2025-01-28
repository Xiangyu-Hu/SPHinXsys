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
 * @file 	kernel_wendland_c2_ck.h
 * @brief 	This is the class for Wendland kernel.
 * @author	Xiangyu Hu
 */

#ifndef KERNEL_WENDLAND_C2_CK_H
#define KERNEL_WENDLAND_C2_CK_H

#include "base_kernel.h"

namespace SPH
{
class KernelWendlandC2CK
{
  public:
    explicit KernelWendlandC2CK(Kernel &kernel)
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
    };

    Real W(const Real &displacement) const
    {
        Real q = displacement * inv_h_;
        return factor_W_1D_ * W_1D(q);
    };

    Real W(const Vec2d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_W_2D_ * W_1D(q);
    };

    Real W(const Vec3d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_W_3D_ * W_1D(q);
    };

    Real W_1D(Real q) const
    {
        return pow(1.0 - 0.5 * q, 4) * (1.0 + 2.0 * q);
    };

    Real dW(const Real &displacement) const
    {
        Real q = displacement * inv_h_;
        return factor_dW_1D_ * dW_1D(q);
    };
    Real dW(const Vec2d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_dW_2D_ * dW_1D(q);
    };
    Real dW(const Vec3d &displacement) const
    {
        Real q = displacement.norm() * inv_h_;
        return factor_dW_3D_ * dW_1D(q);
    };

    Real dW_1D(const Real q) const
    {
        return 0.625 * pow(q - 2.0, 3) * q;
    };

    Vec2d e(const Real &distance, const Vec2d &displacement) const
    {
        return displacement / (distance + TinyReal);
    };
    Vec3d e(const Real &distance, const Vec3d &displacement) const
    {
        return displacement / (distance + TinyReal);
    };

    bool checkIfWithinCutOffRadius(const Vec2d &displacement) const
    {

        return displacement.squaredNorm() < CutOffRadiusSqr();
    };

    bool checkIfWithinCutOffRadius(const Vec3d &displacement) const
    {

        return displacement.squaredNorm() < CutOffRadiusSqr();
    };
    ;

    inline Real CutOffRadius() const { return rc_ref_; };
    inline Real CutOffRadiusSqr() const { return rc_ref_sqr_; };

  private:
    Real inv_h_, rc_ref_, rc_ref_sqr_,
        factor_W_1D_, factor_W_2D_, factor_W_3D_,
        factor_dW_1D_, factor_dW_2D_, factor_dW_3D_;
};
} // namespace SPH
#endif // KERNEL_WENDLAND_C2_CK_H