#ifndef PLASTIC_FUNCTIONS_H
#define PLASTIC_FUNCTIONS_H

#include "data_type.h"
#include <sycl/sycl.hpp>
//=================================================================================================//
namespace SPH
{

//=================================================================================================//
Mat3d upgradeToMat3d_sycl(const Mat2d &input)
{
    Mat3d output = Mat3d::Zero();
    output.block<2, 2>(0, 0) = input;
    return output;
}
//=================================================================================================//
Mat3d upgradeToMat3d_sycl(const Mat3d &input)
{
    return input;
}

//=================================================================================================//
Vec2d degradeToVecd_sycl(const Vec3d &input)
{
    Vec2d output = Vec2d::Zero();
    for (int i = 0; i != Dimensions; i++)
        output[i] = input[i];
    return output;
}
//=================================================================================================//
Mat2d degradeToMatd_sycl(const Mat3d &input)
{
    Mat2d output = Mat2d::Zero();
    for (int i = 0; i != Dimensions; i++)
        for (int j = 0; j != Dimensions; j++)
            output(i, j) = input(i, j);
    return output;
}

Mat3d ConstitutiveRelation_sycl(Mat3d &velocity_gradient, Mat3d &stress_tensor)
{
    // Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    // Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    // Mat3d deviatoric_strain_rate = strain_rate - (1.0 / stress_dimension_) * strain_rate.trace() * Mat3d::Identity();
    // Mat3d stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + K_ * strain_rate.trace() * Mat3d::Identity() + stress_tensor * (spin_rate.transpose()) + spin_rate * stress_tensor;
    // Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
    // Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    // Real f = sqrt(stress_tensor_J2) + alpha_phi_ * stress_tensor.trace() - k_c_;
    // Real lambda_dot_ = 0;
    // Mat3d g = Mat3d::Zero();
    // if (f >= TinyReal)
    // {
    //     Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
    //     // non-associate flow rule
    //     lambda_dot_ = (3.0 * alpha_phi_ * K_ * strain_rate.trace() + (G_ / sqrt(stress_tensor_J2)) * deviatoric_stress_times_strain_rate) / (9.0 * alpha_phi_ * K_ * getDPConstantsA(psi_) + G_);
    //     g = lambda_dot_ * (3.0 * K_ * getDPConstantsA(psi_) * Mat3d::Identity() + G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2)));
    // }
    // Mat3d stress_rate_temp = stress_rate_elastic - g;


    Mat3d stress_rate_temp = Mat3d::Zero();
    return stress_rate_temp;
}


}

#endif // PLASTIC_FUNCTIONS_H
