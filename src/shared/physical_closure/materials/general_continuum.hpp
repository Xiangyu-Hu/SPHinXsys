#ifndef GENERAL_CONTINUUM_HPP
#define GENERAL_CONTINUUM_HPP

#include "general_continuum.h"

namespace SPH
{
//=================================================================================================//
template <typename ExecutionPolicy>
GeneralContinuum::ConstituteKernel::
    ConstituteKernel(const ExecutionPolicy &ex_policy, GeneralContinuum &encloser)
    : E_(encloser.E_), G_(encloser.G_), K_(encloser.K_),
      nu_(encloser.nu_), contact_stiffness_(encloser.contact_stiffness_),
      rho0_(encloser.rho0_) {}
//=================================================================================================//
Real GeneralContinuum::ConstituteKernel::getBulkModulus(Real youngs_modulus, Real poisson_ratio)
{
    return youngs_modulus / 3.0 / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
Real GeneralContinuum::ConstituteKernel::getShearModulus(Real youngs_modulus, Real poisson_ratio)
{
    return 0.5 * youngs_modulus / (1.0 + poisson_ratio);
}
//=================================================================================================//
Real GeneralContinuum::ConstituteKernel::getLambda(Real youngs_modulus, Real poisson_ratio)
{
    return nu_ * youngs_modulus / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
template <typename ExecutionPolicy>
PlasticContinuum::ConstituteKernel::
    ConstituteKernel(const ExecutionPolicy &ex_policy, PlasticContinuum &encloser)
    : GeneralContinuum::ConstituteKernel(ex_policy, encloser),
      c_(encloser.c_), phi_(encloser.phi_),
      psi_(encloser.psi_), alpha_phi_(encloser.alpha_phi_), k_c_(encloser.k_c_) {}
//=================================================================================================//
Real PlasticContinuum::ConstituteKernel::getDPConstantsA(Real friction_angle)
{
    return tan(friction_angle) / sqrt(9.0 + 12.0 * tan(friction_angle) * tan(friction_angle));
};
//=================================================================================================//
Mat3d PlasticContinuum::ConstituteKernel::StressTensorRate(
    UnsignedInt index_i, const Mat3d &velocity_gradient, const Mat3d &stress_tensor)
{
    Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Mat3d deviatoric_strain_rate = strain_rate - (1.0 / stress_dimension_) * strain_rate.trace() * Mat3d::Identity();
    Mat3d stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + K_ * strain_rate.trace() * Mat3d::Identity() +
                                stress_tensor * (spin_rate.transpose()) + spin_rate * stress_tensor;
    Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    Real f = sqrt(stress_tensor_J2) + alpha_phi_ * stress_tensor.trace() - k_c_;
    if (f >= TinyReal)
    {
        Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
        // non_associate flow rule
        Real lambda_dot_ = (3.0 * alpha_phi_ * K_ * strain_rate.trace() +
                            (G_ / (sqrt(stress_tensor_J2) + TinyReal)) * deviatoric_stress_times_strain_rate) /
                           (9.0 * alpha_phi_ * K_ * getDPConstantsA(psi_) + G_);
        return stress_rate_elastic - lambda_dot_ * (3.0 * K_ * getDPConstantsA(psi_) * Mat3d::Identity() +
                                                    G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2 + TinyReal)));
    }
    return stress_rate_elastic;
};
//=================================================================================================//
Mat3d PlasticContinuum::ConstituteKernel::updateStressTensor(
    UnsignedInt index_i, const Mat3d &prev_stress_tensor, const Mat3d &stress_tensor_increment)
{
    return ReturnMapping(index_i, prev_stress_tensor + stress_tensor_increment);
}
//=================================================================================================//
Mat3d PlasticContinuum::ConstituteKernel::ReturnMapping(UnsignedInt index_i, Mat3d stress_tensor)
{
    Real stress_tensor_I1 = stress_tensor.trace();
    if (-alpha_phi_ * stress_tensor_I1 + k_c_ < 0)
        stress_tensor -= (1.0 / stress_dimension_) * (stress_tensor_I1 - k_c_ / alpha_phi_) * Mat3d::Identity();
    stress_tensor_I1 = stress_tensor.trace();
    Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
    volatile Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    if (-alpha_phi_ * stress_tensor_I1 + k_c_ < sqrt(stress_tensor_J2))
    {
        Real r = (-alpha_phi_ * stress_tensor_I1 + k_c_) / (sqrt(stress_tensor_J2) + TinyReal);
        stress_tensor = r * deviatoric_stress_tensor + (1.0 / stress_dimension_) * stress_tensor_I1 * Mat3d::Identity();
    }
    return stress_tensor;
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_CONTINUUM_HPP