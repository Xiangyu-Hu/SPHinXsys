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
 * @file 	base_kernel.h
 * @brief 	This is the base classes of kernel functions.  Implementation will be
 *			implemented in derived classes. The kernal function define the relevance
 * 			between two neighboring particles. Basically, the further the two
 *			particles, the less relevance they have.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_KERNELS_H
#define BASE_KERNELS_H

#include "base_data_package.h"

#include <functional>
#include <string>

using namespace std::placeholders;

namespace SPH
{
/**
 * @class 	Kernel
 * @brief 	Abstract base class of a general SPH kernel function which
 * 			is a smoothed Dirac delta function,
 * 			a kernel function is radial symmetric, and has a scaling factor.
 * 			Based on difference data type in 2d or 3d buildings,
 * 			the kernel is defined for 2 and 3 dimensions.
 * 			The kernel gives value one at the origin.
 * 			The naming of kernel function follows the stand SPH literature.
 * 			Currently, only constant smoothing length is applied.
 * 			Basically, one can assign different kernel for different particle interactions.
 */
class Kernel
{
  protected:
    std::string kernel_name_;
    Real h_, inv_h_;   /**< reference smoothing length and its inverse **/
    Real kernel_size_; /**<kernel_size_ *  h_ gives the zero kernel value */
    Real truncation_;  /**< to obtain cut off radius */
    Real rc_ref_;      /** reference cut off radius, beyond this kernel value is neglected. **/
    Real rc_ref_sqr_;  /** reference cut off radius square **/
    /** Normalization factors for the kernel function  **/
    Real factor_W_1D_, factor_W_2D_, factor_W_3D_;
    /** Auxiliary factors for the derivative of kernel function  **/
    Real factor_dW_1D_, factor_dW_2D_, factor_dW_3D_;
    /** Auxiliary factors for the second order derivative of kernel function  **/
    Real factor_d2W_1D_, factor_d2W_2D_, factor_d2W_3D_;

    void setDerivativeParameters();

  public:
    /** empty initialization in constructor, initialization will be carried out later. */
    explicit Kernel(Real h, Real kernel_size, Real rc_ref, const std::string &name);
    virtual ~Kernel(){};

    std::string Name() const { return kernel_name_; };
    void resetSmoothingLength(Real h);
    Real SmoothingLength() const { return h_; };
    /**< non-dimensional size of the kernel, generally 2.0 **/
    Real KernelSize() const { return kernel_size_; };
    Real Truncation() const { return truncation_; };
    Real CutOffRadius() const { return rc_ref_; };
    Real CutOffRadiusSqr() const { return rc_ref_sqr_; };
    Real FactorW1D() const { return factor_W_1D_; };
    Real FactorW2D() const { return factor_W_2D_; };
    Real FactorW3D() const { return factor_W_3D_; };

    /**
     * Calculates the kernel value for the given displacement of two particles
     * r_ij pointing from particle j to particle i
     */
    Real W(const Real &r_ij, const Real &displacement) const;
    Real W(const Real &r_ij, const Vec2d &displacement) const;
    Real W(const Real &r_ij, const Vec3d &displacement) const;

    /** this value could be use to calculate the value of W
     * they are realized in specific kernel implementations
     */
    virtual Real W_1D(const Real q) const = 0;
    virtual Real W_2D(const Real q) const = 0;
    virtual Real W_3D(const Real q) const = 0;

    /** Calculates the kernel value at the origin **/
    Real W0(const Real &point_i) const { return factor_W_1D_; };
    Real W0(const Vec2d &point_i) const { return factor_W_2D_; };
    Real W0(const Vec3d &point_i) const { return factor_W_3D_; };

    /** Calculates the kernel derivation for
     * the given distance of two particles
     */
    Real dW(const Real &r_ij, const Real &displacement) const;
    Real dW(const Real &r_ij, const Vec2d &displacement) const;
    Real dW(const Real &r_ij, const Vec3d &displacement) const;

    /** this value could be use to calculate the value of dW
     * they are realized in specific kernel implementations
     */
    virtual Real dW_1D(const Real q) const = 0;
    virtual Real dW_2D(const Real q) const = 0;
    virtual Real dW_3D(const Real q) const = 0;

    /** Calculates the kernel second order derivation for
     * the given distance of two particles
     */
    Real d2W(const Real &r_ij, const Real &displacement) const;
    Real d2W(const Real &r_ij, const Vec2d &displacement) const;
    Real d2W(const Real &r_ij, const Vec3d &displacement) const;

    /** this value could be use to calculate the value of d2W
     * they are realized in specific kernel implementations
     */
    virtual Real d2W_1D(const Real q) const = 0;
    virtual Real d2W_2D(const Real q) const = 0;
    virtual Real d2W_3D(const Real q) const = 0;

    //----------------------------------------------------------------------
    //		Below are for variable smoothing length.
    //		Note that we input the ratio between the reference smoothing length
    //		to the variable smoothing length.
    //----------------------------------------------------------------------
  protected:
    /** Functor for variable smoothing length. */
    typedef std::function<Real(const Real &)> FactorFunctor;
    FactorFunctor h_factor_W_1D_, h_factor_W_2D_, h_factor_W_3D_;
    FactorFunctor h_factor_dW_1D_, h_factor_dW_2D_, h_factor_dW_3D_;
    FactorFunctor h_factor_d2W_1D_, h_factor_d2W_2D_, h_factor_d2W_3D_;

    Real factorW1D(const Real &h_ratio) const { return h_ratio; };
    Real factorW2D(const Real &h_ratio) const { return h_ratio * h_ratio; };
    Real factorW3D(const Real &h_ratio) const { return h_ratio * h_ratio * h_ratio; };
    Real factordW1D(const Real &h_ratio) const { return factorW1D(h_ratio) * h_ratio; };
    Real factordW2D(const Real &h_ratio) const { return factorW2D(h_ratio) * h_ratio; };
    Real factordW3D(const Real &h_ratio) const { return factorW3D(h_ratio) * h_ratio; };
    Real factord2W1D(const Real &h_ratio) const { return factordW1D(h_ratio) * h_ratio; };
    Real factord2W2D(const Real &h_ratio) const { return factordW2D(h_ratio) * h_ratio; };
    Real factord2W3D(const Real &h_ratio) const { return factordW3D(h_ratio) * h_ratio; };

  public:
    Real CutOffRadius(Real h_ratio) const { return rc_ref_ / h_ratio; };
    Real CutOffRadiusSqr(Real h_ratio) const { return rc_ref_sqr_ / (h_ratio * h_ratio); };

    Real W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const;
    Real W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const;
    Real W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const;

    /** Calculates the kernel value at the origin **/
    Real W0(const Real &h_ratio, const Real &point_i) const;
    Real W0(const Real &h_ratio, const Vec2d &point_i) const;
    Real W0(const Real &h_ratio, const Vec3d &point_i) const;

    /** Calculates the kernel derivation for the given distance of two particles **/
    Real dW(const Real &h_ratio, const Real &r_ij, const Real &displacement) const;
    Real dW(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const;
    Real dW(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const;

    /** Calculates the kernel second order derivation for the given distance of two particles **/
    Real d2W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const;
    Real d2W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const;
    Real d2W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const;
    //----------------------------------------------------------------------
    //		Below are for reduced kernels.
    //----------------------------------------------------------------------
  public:
    void reduceOnce();  /** reduce for thin structures or films */
    void reduceTwice(); /** reduce for linear structures or filaments */
};
} // namespace SPH
#endif // BASE_KERNELS_H