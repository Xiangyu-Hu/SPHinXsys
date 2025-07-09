#ifndef GENERAL_CONTINUUM_HPP
#define GENERAL_CONTINUUM_HPP

#include "general_continuum.h"

namespace SPH
{
Real GeneralContinuum::GeneralContinuumKernel::getBulkModulus(Real youngs_modulus, Real poisson_ratio)
{
    return youngs_modulus / 3.0 / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
Real GeneralContinuum::GeneralContinuumKernel::getShearModulus(Real youngs_modulus, Real poisson_ratio)
{
    return 0.5 * youngs_modulus / (1.0 + poisson_ratio);
}
//=================================================================================================//
Real GeneralContinuum::GeneralContinuumKernel::getLambda(Real youngs_modulus, Real poisson_ratio)
{
    return nu_ * youngs_modulus / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
}
//=================================================================================================//
//=================================================================================================//
Real PlasticContinuum::PlasticKernel::getDPConstantsA(Real friction_angle)
{
    return tan(friction_angle) / sqrt(9.0 + 12.0 * tan(friction_angle) * tan(friction_angle));
};
//=================================================================================================//
Mat3d PlasticContinuum::PlasticKernel::ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor)
{
    Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Mat3d deviatoric_strain_rate = strain_rate - (1.0 / stress_dimension_) * strain_rate.trace() * Mat3d::Identity();
    Mat3d stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + K_ * strain_rate.trace() * Mat3d::Identity() + stress_tensor * (spin_rate.transpose()) + spin_rate * stress_tensor;
    Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    Real f = sqrt(stress_tensor_J2) + alpha_phi_ * stress_tensor.trace() - k_c_;
    Real lambda_dot_ = 0;
    Mat3d g = Mat3d::Zero();
    if (f >= TinyReal)
    {
        Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
        // non-associate flow rule
        lambda_dot_ = (3.0 * alpha_phi_ * K_ * strain_rate.trace() + (G_ / (sqrt(stress_tensor_J2)+TinyReal)) * deviatoric_stress_times_strain_rate) / (9.0 * alpha_phi_ * K_ * getDPConstantsA(psi_) + G_);
        g = lambda_dot_ * (3.0 * K_ * getDPConstantsA(psi_) * Mat3d::Identity() + G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2+ TinyReal)));
    }
    Mat3d stress_rate_temp = stress_rate_elastic - g;
    return stress_rate_temp;
}; 
//=================================================================================================//
Mat3d PlasticContinuum::PlasticKernel::ReturnMapping(Mat3d &stress_tensor)
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
Real PlasticContinuum::PlasticKernel::getFrictionVelocity(Real uz, Real z)
{
    const Real nu = 1.0e-6;        // Kinematic viscosity (m^2/s)
    const Real d50 = d_s_;         // Median grain diameter (m)
    const Real kappa = 0.41;       // von Kármán constant
    const Real ks = 2.0 * d50;     // Equivalent roughness height (m)
    const Real tol = 0.01;         // Relative error tolerance

    // Initial estimate of shear velocity u_star assuming laminar sublayer
    Real u_star_1 = std::sqrt(uz * nu / z);
    Real delta_v = 11.6 * nu / u_star_1;

    // Check smooth flow condition
    if ((ks * u_star_1 / nu < 5.0) && (z < delta_v)) {
        return u_star_1;
    }

    // Iterative solution for u_star in rough/turbulent flow
    Real u_star_old = u_star_1;
    Real u_star_new = 0.0;
    int iter = 0;

    while (true) {
        // Compute roughness length z0 based on u_star
        Real ksu_nu = ks * u_star_old / nu;
        Real z0;

        if (ksu_nu < 5.0) {
            z0 = 0.11 * nu / u_star_old;
        } else if (ksu_nu > 70.0) {
            z0 = 0.033 * ks;
        } else {
            z0 = 0.11 * nu / u_star_old + 0.033 * ks;
        }

        // Update u_star using the logarithmic velocity profile
        u_star_new = (kappa * uz) / std::log(z / z0);

        // Check convergence
        Real rel_err = std::abs(u_star_new - u_star_old) / u_star_old;
        if (rel_err < tol || iter > 100) break;

        u_star_old = u_star_new;
        ++iter;
    }

    return u_star_new;
}
//=================================================================================================//
Real PlasticContinuum::PlasticKernel::calculateThetaCr(Real u_star)
{
    const Real rho_w = 1000;  // 
    const Real mu_w = 0.001;  // 

    Real Re = rho_w * u_star * d_s_ / mu_w;
    if (Re <= 500 && Re > 0.0) 
    {
        return 0.010595 * log(Re) + 0.110476 / Re + 0.0027197;
    } else 
    {
        return 0.068;
    } 
}
//=================================================================================================//
Real PlasticContinuum::PlasticKernel::ThetaToFrictionVelcoty(Real theta_cr)
{
    const Real rho_w = 1000.0;        // Water density (kg/m³)
    const Real gravity = 9.8;         // Gravitational acceleration (m/s²)
    const Real Sn = rho0_ / rho_w;    // Relative density = particle density / water density

    // Compute critical shear velocity using:
    // u*_cr = sqrt[ θ_cr * ( (ρ_s - ρ_w) * g * d ) / ρ_w ]
    Real friction_velocity = sqrt(theta_cr * ((Sn - 1.0) * gravity * d_s_) / 1.0);

    return friction_velocity;
}
//=================================================================================================//
}// namespace SPH
#endif //GENERAL_CONTINUUM_HPP