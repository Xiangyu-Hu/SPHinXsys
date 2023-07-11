/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2023 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	anisotropic_kernel.h
 * @brief 	This is the base classes of anisotropic kernel functions.  Implementation will be
 *			implemented in derived classes. The kernel function define the relevance
 * 			between two neighboring particles. Basically, the further the two
 *			particles, the less relevance they have.
 * @author	 Zhentong Wang and Xiaojing Tang
 * @version	0.1
 * @version  0.3.0
 */

#ifndef ANISOTROPIC_KERNEL_H
#define ANISOTROPIC_KERNEL_H

#include "base_kernel_includes_nonisotropic.h"
#include "kernel_wenland_c2_anisotropic.h"

namespace SPH
{

/**
 * @class AnisotropicKernel
 * @brief Abstract base class of a general anisotropic  SPH kernel function which
 * is a smoothed Dirac delta function.
 * Based on difference data type in 2d or 3d buildings,
 * the kernel is defined for 2 and 3 dimensions.
 * The kernel gives value one at the origin.
 * The naming of kernel function follows the stand SPH literature.
 * Basically, one can assign different kernel for different particle interactions.
 */

template <class KernelType>
class AnisotropicKernel : public KernelType
{
  protected:
    Mat2d transformed_tensor_2d_;
    Mat3d transformed_tensor_3d_;

  public:
    explicit AnisotropicKernel(Real h, Vec2d kernel_vector = Vec2d(1.0, 1.0), Vec2d transform_vector = Vec2d(0.0, 0.0))
        : KernelType(h), transformed_tensor_2d_(1.0 / h * Mat2d::Identity()), transformed_tensor_3d_(1.0 / h * Mat3d::Identity())
    {
        this->kernel_name_ = "AnisotropicKernel_" + this->kernel_name_;
        transformed_tensor_2d_ = getCoordinateTransformationTensorG(kernel_vector, transform_vector);
        getFactors();
    };

    explicit AnisotropicKernel(Real h, Vec3d kernel_vector = Vec3d(1.0, 1.0, 1.0), Vec3d transform_vector = Vec3d::Zero())
        : KernelType(h), transformed_tensor_2d_(1.0 / h * Mat2d::Identity()), transformed_tensor_3d_(1.0 / h * Mat3d::Identity())
    {
        this->kernel_name_ = "AnisotropicKernel_" + this->kernel_name_;
        transformed_tensor_3d_ = getCoordinateTransformationTensorG(kernel_vector, transform_vector);
        getFactors();
    };

    void getFactors()
    {
        this->factor_W_1D_ = this->h_ * 1.0 / this->h_ * this->FactorW1D();
        this->factor_W_2D_ = this->h_ * this->h_ * transformed_tensor_2d_.determinant() * this->FactorW2D();
        this->factor_W_3D_ = this->h_ * this->h_ * this->h_ * transformed_tensor_3d_.determinant() * this->FactorW3D();

        this->factor_dW_1D_ = this->factor_W_1D_;
        this->factor_dW_2D_ = this->factor_W_2D_;
        this->factor_dW_3D_ = this->factor_W_3D_;

        this->factor_d2W_1D_ = this->factor_W_1D_;
        this->factor_d2W_2D_ = this->factor_W_2D_;
        this->factor_d2W_3D_ = this->factor_W_3D_;
    };

    Mat2d getCoordinateTransformationTensorG(Vec2d kernel_vector, Vec2d transform_vector);
    Mat3d getCoordinateTransformationTensorG(Vec3d kernel_vector, Vec3d transform_vector);

    virtual Vec2d e(const Real &distance, const Vec2d &displacement) const override;
    virtual Vec3d e(const Real &distance, const Vec3d &displacement) const override;

    virtual bool checkIfWithinCutOffRadius(Vec2d displacement) override;
    virtual bool checkIfWithinCutOffRadius(Vec3d displacement) override;

    /** Calculates the kernel value at the origin **/
    virtual Real W0(const Real &point_i) const override { return this->factor_W_1D_; };
    virtual Real W0(const Vec2d &point_i) const override { return this->factor_W_2D_; };
    virtual Real W0(const Vec3d &point_i) const override { return this->factor_W_3D_; };

    /** Calculates the kernel value for the given displacement of two particles
     * r_ij pointing from particle j to particle i
     */
    virtual Real W(const Real &r_ij, const Real &displacement) const override;
    virtual Real W(const Real &r_ij, const Vec2d &displacement) const override;
    virtual Real W(const Real &r_ij, const Vec3d &displacement) const override;

    /** Calculates the kernel derivation for
     * the given distance of two particles
     */
    virtual Real dW(const Real &r_ij, const Real &displacement) const override;
    virtual Real dW(const Real &r_ij, const Vec2d &displacement) const override;
    virtual Real dW(const Real &r_ij, const Vec3d &displacement) const override;

    /** Calculates the kernel second order derivation for
     * the given distance of two particles
     */
    virtual Real d2W(const Real &r_ij, const Real &displacement) const override;
    virtual Real d2W(const Real &r_ij, const Vec2d &displacement) const override;
    virtual Real d2W(const Real &r_ij, const Vec3d &displacement) const override;

    virtual ~AnisotropicKernel(){};
};

} // namespace SPH
#endif // ANISOTROPIC_KERNEL_H
