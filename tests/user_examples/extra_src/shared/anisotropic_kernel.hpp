/**
 * @file 	anisotropic_kernel.hpp
 * @brief
 * @author
 */
#ifndef ANISOTROPIC_KERNEL_HPP
#define ANISOTROPIC_KERNEL_HPP

#include "anisotropic_kernel.h"

namespace SPH
{
//=========================================================================================//
template <class KernelType>
Mat2d AnisotropicKernel<KernelType>::getCoordinateTransformationTensorG(Vec2d kernel_vector, Vec2d transform_angle)
{
    Mat2d scaling_tensor = Mat2d({{ Real(1.0) / (this->h_ * kernel_vector[0]), 0.0},
                                  {0.0, Real(1.0) / (this->h_ * kernel_vector[1])}}
                                  );
    Mat2d G_kernel_coordinate = scaling_tensor;
    return G_kernel_coordinate;
}
//=========================================================================================//
template <class KernelType>
Mat3d AnisotropicKernel<KernelType>::getCoordinateTransformationTensorG(Vec3d kernel_vector, Vec3d transform_vector)
{
    Mat3d G_kernel_coordinate = Mat3d({{ Real(1.0) / (this->h_ * kernel_vector[0]), 0.0, 0.0},
                                       {0.0, Real(1.0) / (this->h_ * kernel_vector[1]), 0.0},
                                       {0.0, 0.0, Real(1.0) / (this->h_ * kernel_vector[2])}}
                                     );
    return G_kernel_coordinate;
}
//=========================================================================================//
template <class KernelType>
Vec2d AnisotropicKernel<KernelType>::e(const Real &distance, const Vec2d &displacement) const
{
    Vec2d transformed_displacement_ = this->h_ * transformed_tensor_2d_ * displacement;
    return transformed_tensor_2d_ * transformed_displacement_ / ((transformed_displacement_).norm() + TinyReal);
}
//=========================================================================================//
template <class KernelType>
Vec3d AnisotropicKernel<KernelType>::e(const Real &r_ij, const Vec3d &displacement) const
{
    Vec3d transformed_displacement_ = this->h_ * transformed_tensor_3d_ * displacement;
    return transformed_tensor_3d_ * transformed_displacement_ / ((transformed_displacement_).norm() + TinyReal);
}
//=========================================================================================//
template <class KernelType>
bool AnisotropicKernel<KernelType>::checkIfWithinCutOffRadius(Vec2d displacement)
{
    Vec2d transformed_displacement = this->h_ * transformed_tensor_2d_ * displacement;
    Real distance_metric = transformed_displacement.squaredNorm();
    if (distance_metric < this->CutOffRadiusSqr())
        return true;
    else
        return false;
}
//=========================================================================================//
template <class KernelType>
bool AnisotropicKernel<KernelType>::checkIfWithinCutOffRadius(Vec3d displacement)
{
    Vec3d transformed_displacement = this->h_ * transformed_tensor_3d_ * displacement;
    Real distance_metric = transformed_displacement.squaredNorm();
    if (distance_metric < this->CutOffRadiusSqr())
        return true;
    else
        return false;
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::W(const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij / this->h_;
    return this->factor_W_1D_ * this->W_1D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::W(const Real &r_ij, const Vec2d &displacement) const
{
    Vec2d transformed_displacement = transformed_tensor_2d_ * displacement;
    Real q = transformed_displacement.norm();
    return this->factor_W_2D_ * this->W_2D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::W(const Real &r_ij, const Vec3d &displacement) const
{
    Vec3d transformed_displacement = transformed_tensor_3d_ * displacement;
    Real q = transformed_displacement.norm();
    return this->factor_W_3D_ * this->W_3D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::dW(const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij / this->h_;
    return this->factor_dW_1D_ * this->dW_1D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::dW(const Real &r_ij, const Vec2d &displacement) const
{
    Vec2d transformed_displacement = transformed_tensor_2d_ * displacement;
    Real q = transformed_displacement.norm();
    return this->factor_dW_2D_ * this->dW_2D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::dW(const Real &r_ij, const Vec3d &displacement) const
{
    Vec3d transformed_displacement = transformed_tensor_3d_ * displacement;
    Real q = transformed_displacement.norm();
    return this->factor_dW_3D_ * this->dW_3D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::d2W(const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij / this->h_;
    return this->factor_d2W_1D_ * this->d2W_1D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::d2W(const Real &r_ij, const Vec2d &displacement) const
{
    Vec2d transformed_displacement = transformed_tensor_2d_ * displacement;
    Real q = transformed_displacement.norm();
    Real derivate_parameter2D_ = (transformed_tensor_2d_ * transformed_displacement).norm() / (q + TinyReal);
    return this->factor_d2W_2D_ * derivate_parameter2D_ * derivate_parameter2D_ * this->d2W_2D(q);
}
//=========================================================================================//
template <class KernelType>
Real AnisotropicKernel<KernelType>::d2W(const Real &r_ij, const Vec3d &displacement) const
{
    Vec3d transformed_displacement = transformed_tensor_3d_ * displacement;
    Real q = transformed_displacement.norm();
    Real derivate_parameter3D_ = (transformed_tensor_3d_ * transformed_displacement).norm() / (q + TinyReal);
    return this->factor_d2W_3D_ * derivate_parameter3D_ * derivate_parameter3D_ * this->d2W_3D(q);
}
//=========================================================================================//
} // namespace SPH
#endif
