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
Mat3d PlasticContinuum::ConstituteKernel::
    updateStressTensor(UnsignedInt index_i, const Mat3d &try_stress_tensor)
{
    return ReturnMapping(index_i, try_stress_tensor);
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
template <typename ExecutionPolicy>
J2Plasticity::ConstituteKernel::
    ConstituteKernel(const ExecutionPolicy &ex_policy, J2Plasticity &encloser)
    : GeneralContinuum::ConstituteKernel(ex_policy, encloser),
      yield_stress_(encloser.yield_stress_),
      hardening_modulus_(encloser.hardening_modulus_),
      sqrt_2_over_3_(encloser.sqrt_2_over_3_),
      hardening_factor_(encloser.dv_hardening_factor_->DelegatedData(ex_policy)) {}
//=================================================================================================//
Matd J2Plasticity::ConstituteKernel::ShearStressRate(
    UnsignedInt index_i, const Matd &velocity_gradient, const Matd &shear_stress)
{
    Matd strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Matd deviatoric_strain_rate = strain_rate - (1.0 / (Real)Dimensions) * strain_rate.trace() * Matd::Identity();
    Matd spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Matd shear_stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + shear_stress * (spin_rate.transpose()) + spin_rate * shear_stress;
    Real stress_tensor_J2 = 0.5 * (shear_stress.cwiseProduct(shear_stress.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor_[index_i] + yield_stress_);
    if (f > TinyReal)
    {
        Real deviatoric_stress_times_strain_rate = (shear_stress.cwiseProduct(strain_rate)).sum();
        Real lambda_dot_ = deviatoric_stress_times_strain_rate / (sqrt(2.0 * stress_tensor_J2) * (1.0 + hardening_modulus_ / (3.0 * G_)));
        return shear_stress_rate_elastic - lambda_dot_ * (sqrt(2.0) * G_ * shear_stress / (sqrt(stress_tensor_J2)));
    }
    return shear_stress_rate_elastic;
};
//=================================================================================================//
Matd J2Plasticity::ConstituteKernel::updateShearStress(
    UnsignedInt index_i, const Matd &try_shear_stress)
{
    Real hardening_factor_increment = HardeningFactorRate(try_shear_stress, hardening_factor_[index_i]);
    hardening_factor_[index_i] += sqrt(2.0 / 3.0) * hardening_factor_increment;
    return ReturnMapping(index_i, try_shear_stress);
}
//=================================================================================================//
Real J2Plasticity::ConstituteKernel::ScalePenaltyForce(UnsignedInt index_i, const Matd &try_shear_stress)
{
    Real stress_tensor_J2 = 0.5 * (try_shear_stress.cwiseProduct(try_shear_stress.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor_[index_i] + yield_stress_);
    return (f > TinyReal)
               ? (sqrt_2_over_3_ *
                  (hardening_modulus_ * hardening_factor_[index_i] + yield_stress_)) /
                     (sqrt(2.0 * stress_tensor_J2) + TinyReal)
               : 1.0;
}
//=================================================================================================//
Matd J2Plasticity::ConstituteKernel::ReturnMapping(UnsignedInt index_i, Matd try_shear_stress)
{
    return ScalePenaltyForce(index_i, try_shear_stress) * try_shear_stress;
}
//=================================================================================================//
Real J2Plasticity::ConstituteKernel::HardeningFactorRate(
    const Matd &shear_stress, Real &hardening_factor)
{
    Real stress_tensor_J2 = 0.5 * (shear_stress.cwiseProduct(shear_stress.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor + yield_stress_);
    return (f > TinyReal) ? 0.5 * f / (G_ + hardening_modulus_ / 3.0) : 0.0;
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_CONTINUUM_HPP