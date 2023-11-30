#include "general_continuum.h"

namespace SPH
{
//==============================GeneralContinuum===============================================//
Real GeneralContinuum::getBulkModulus(Real youngs_modulus, Real poisson_ratio)
{
    return youngs_modulus / 3.0 / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
Real GeneralContinuum::getShearModulus(Real youngs_modulus, Real poisson_ratio)
{
    return 0.5 * youngs_modulus / (1.0 + poisson_ratio);
}
//=================================================================================================//
Real GeneralContinuum::getLambda(Real youngs_modulus, Real poisson_ratio)
{
    return nu_ * youngs_modulus / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
Matd GeneralContinuum::ConstitutiveRelationShearStress(Matd &velocity_gradient, Matd &shear_stress)
{
    Matd strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Matd spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Matd deviatoric_strain_rate = strain_rate - (1.0 / (Real)Dimensions) * strain_rate.trace() * Matd::Identity();
    Matd stress_rate = 2.0 * G_ * deviatoric_strain_rate + shear_stress * (spin_rate.transpose()) + spin_rate * shear_stress;
    return stress_rate;
}
//=============================================================================================//
//==============================PlasticContinuum===============================================//
//=============================================================================================//

Real PlasticContinuum::getDPConstantsA(Real friction_angle)
{
    return tan(friction_angle) / sqrt(9 + 12 * tan(friction_angle) * tan(friction_angle));
}

Real PlasticContinuum::getDPConstantsK(Real cohesion, Real friction_angle)
{
    return 3 * cohesion / sqrt(9 + 12 * tan(friction_angle) * tan(friction_angle));
}

Mat3d PlasticContinuum::ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor)
{
    Real dim = 3;
    Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Mat3d deviatoric_strain_rate = strain_rate - (1.0 / dim) * strain_rate.trace() * Mat3d::Identity();
    // consider the elastic part
    Mat3d stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + K_ * strain_rate.trace() * Mat3d::Identity() + stress_tensor * (spin_rate.transpose()) + spin_rate * stress_tensor;
    // consider the plastic part
    Real stress_tensor_I1 = stress_tensor.trace(); // first invariant of stress
    Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / dim) * stress_tensor.trace() * Mat3d::Identity();
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    Real f = sqrt(stress_tensor_J2) + alpha_phi_ * stress_tensor_I1 - k_c_;
    Real lambda_dot_ = 0;
    Mat3d g = Mat3d::Zero();
    if ((f >= TinyReal) && (stress_tensor_J2 > TinyReal))
    {
        Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
        // non-associate flow rule
        lambda_dot_ = (3 * alpha_phi_ * K_ * strain_rate.trace() + (G_ / sqrt(stress_tensor_J2)) * deviatoric_stress_times_strain_rate) / (9 * alpha_phi_ * K_ * getDPConstantsA(psi_) + G_);
        g = lambda_dot_ * (3 * K_ * getDPConstantsA(psi_) * Mat3d::Identity() + G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2)));
    }
    Mat3d stress_rate_temp = stress_rate_elastic - g;
    return stress_rate_temp;
}

Mat3d PlasticContinuum::ReturnMapping(Mat3d &stress_tensor)
{
    Real dim = 3;

    Real stress_tensor_I1 = stress_tensor.trace();
    if (-alpha_phi_ * stress_tensor_I1 + k_c_ < 0)
    {
        stress_tensor -= (1.0 / dim) * (stress_tensor_I1 - k_c_ / alpha_phi_) * Mat3d::Identity();
    }
    stress_tensor_I1 = stress_tensor.trace();
    Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / dim) * stress_tensor.trace() * Mat3d::Identity();
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    if (-alpha_phi_ * stress_tensor_I1 + k_c_ < sqrt(stress_tensor_J2))
    {
        Real r = (-alpha_phi_ * stress_tensor_I1 + k_c_) / (sqrt(stress_tensor_J2) + TinyReal);
        stress_tensor = r * deviatoric_stress_tensor + (1.0 / dim) * stress_tensor_I1 * Mat3d::Identity();
    }
    return stress_tensor;
}

//=================================================================================================//
//=====================================J2Plasticity================================================//
//=================================================================================================//
Mat3d J2Plasticity::ConstitutiveRelationShearStress(const Mat3d &velocity_gradient, const Mat3d &shear_stress)
{
    Real dim = 3.0;
    Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    //    Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Mat3d deviatoric_strain_rate = strain_rate - (1.0 / dim) * strain_rate.trace() * Mat3d::Identity();
    // consider the elastic part
    Mat3d shear_stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate;
    // consider the plastic part
    Mat3d deviatoric_stress_tensor = shear_stress;
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt(2.0 / 3.0) * yield_stress_;
    Real lambda_dot_ = 0;
    Mat3d g = Mat3d::Zero();
    if (f > TinyReal)
    {
        Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
        lambda_dot_ = deviatoric_stress_times_strain_rate / sqrt(2.0 * stress_tensor_J2);
        g = lambda_dot_ * (sqrt(2.0) * G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2)));
    }
    return shear_stress_rate_elastic;//    -g;
}
//=================================================================================================//
Mat3d J2Plasticity::ConstitutiveRelationShearStress(Mat3d &velocity_gradient, Mat3d &shear_stress, Real hardening_parameter)
{
    Real dim = 3.0;
    Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    //        Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Mat3d deviatoric_strain_rate = strain_rate - (1.0 / dim) * strain_rate.trace() * Mat3d::Identity();
    // consider the elastic part
    Mat3d shear_stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate;
    // consider the plastic part
    Mat3d deviatoric_stress_tensor = shear_stress;
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    // Real f = sqrt(2.0 * stress_tensor_J2) - sqrt(2.0 / 3.0) * yield_stress_;
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt(2.0 / 3.0) * (hardening_modulus_ * hardening_parameter + yield_stress_);
    Real lambda_dot_ = 0;
    Mat3d g = Mat3d::Zero();
    if (f > TinyReal)
    {
        Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
        lambda_dot_ = deviatoric_stress_times_strain_rate / sqrt(2.0 * stress_tensor_J2);
        g = lambda_dot_ * (sqrt(2.0) * G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2)));
    }
    return shear_stress_rate_elastic - g;
}
//=================================================================================================//
Mat3d J2Plasticity::ReturnMappingShearStress(const Mat3d &shear_stress)
{
    Mat3d deviatoric_stress_tensor = shear_stress;
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    Real r = 0;
    if (sqrt(2.0 * stress_tensor_J2) - sqrt(2.0 / 3.0) * yield_stress_ > TinyReal)
    {
        r = (sqrt(2.0 / 3.0) * yield_stress_) / (sqrt(2.0 * stress_tensor_J2) + TinyReal);
        deviatoric_stress_tensor *= r;
    }
    return deviatoric_stress_tensor;
}
//=================================================================================================//
Mat3d J2Plasticity::ReturnMappingShearStress(Mat3d &shear_stress, Real hardening_parameter)
{
    Mat3d deviatoric_stress_tensor = shear_stress;
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    Real r = 0;
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt(2.0 / 3.0) * (hardening_modulus_ * hardening_parameter + yield_stress_);
    if (f > TinyReal)
    {

        r = (sqrt(2.0 / 3.0) * (hardening_modulus_ * hardening_parameter + yield_stress_)) / (sqrt(2.0 * stress_tensor_J2) + TinyReal);
        shear_stress = r * deviatoric_stress_tensor;
    }
    return shear_stress;
}
//=================================================================================================//
int J2Plasticity::PlasticIndicator(const Mat3d &shear_stress)
{
    int indicator = 0;
    Mat3d deviatoric_stress_tensor = shear_stress;
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    if (sqrt(2.0 * stress_tensor_J2) - sqrt(2.0 / 3.0) * yield_stress_ > TinyReal)
    {
        indicator = 1;
    }
    return indicator;
}
} // namespace SPH