#include "general_continuum.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
GeneralContinuum::GeneralContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio)
    : WeaklyCompressibleFluid(rho0, c0), SolidContact(rho0, rho0 * c0 * c0, 0.0),
      E_(0.0), G_(0.0), K_(0.0), nu_(0.0)
{
    material_type_name_ = "GeneralContinuum";
    E_ = youngs_modulus;
    nu_ = poisson_ratio;
    G_ = getShearModulus(youngs_modulus, poisson_ratio);
    K_ = getBulkModulus(youngs_modulus, poisson_ratio);
}
//=================================================================================================//
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
//=================================================================================================//
PlasticContinuum::PlasticContinuum(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio,
                                   Real friction_angle, Real cohesion, Real dilatancy)
    : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
      c_(cohesion), phi_(friction_angle), psi_(dilatancy), alpha_phi_(0.0), k_c_(0.0)
{
    material_type_name_ = "PlasticContinuum";
    alpha_phi_ = getDPConstantsA(friction_angle);
    k_c_ = getDPConstantsK(cohesion, friction_angle);
}
//=================================================================================================//
Real PlasticContinuum::getDPConstantsA(Real friction_angle)
{
    return tan(friction_angle) / sqrt(9.0 + 12.0 * tan(friction_angle) * tan(friction_angle));
}
//=================================================================================================//
Real PlasticContinuum::getDPConstantsK(Real cohesion, Real friction_angle)
{
    return 3.0 * cohesion / sqrt(9.0 + 12.0 * tan(friction_angle) * tan(friction_angle));
}
//=================================================================================================//
Mat3d PlasticContinuum::ConstitutiveRelation(Mat3d &velocity_gradient, Mat3d &stress_tensor)
{
    Mat3d strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Mat3d spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Mat3d deviatoric_strain_rate = strain_rate - (1.0 / stress_dimension_) * strain_rate.trace() * Mat3d::Identity();
    Mat3d stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + K_ * strain_rate.trace() * Mat3d::Identity() +
                                stress_tensor * (spin_rate.transpose()) + spin_rate * stress_tensor;
    Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    Real f = sqrt(stress_tensor_J2) + alpha_phi_ * stress_tensor.trace() - k_c_;
    Real lambda_dot_ = 0;
    Mat3d g = Mat3d::Zero();
    if (f >= TinyReal)
    {
        Real deviatoric_stress_times_strain_rate = (deviatoric_stress_tensor.cwiseProduct(strain_rate)).sum();
        // non-associate flow rule
        lambda_dot_ = (3.0 * alpha_phi_ * K_ * strain_rate.trace() + (G_ / sqrt(stress_tensor_J2)) * deviatoric_stress_times_strain_rate) /
                      (9.0 * alpha_phi_ * K_ * getDPConstantsA(psi_) + G_);
        g = lambda_dot_ * (3.0 * K_ * getDPConstantsA(psi_) * Mat3d::Identity() + G_ * deviatoric_stress_tensor / (sqrt(stress_tensor_J2)));
    }
    Mat3d stress_rate_temp = stress_rate_elastic - g;
    return stress_rate_temp;
}
//=================================================================================================//
Mat3d PlasticContinuum::ReturnMapping(Mat3d &stress_tensor)
{
    Real stress_tensor_I1 = stress_tensor.trace();
    if (-alpha_phi_ * stress_tensor_I1 + k_c_ < 0)
        stress_tensor -= (1.0 / stress_dimension_) * (stress_tensor_I1 - k_c_ / alpha_phi_) * Mat3d::Identity();
    stress_tensor_I1 = stress_tensor.trace();
    Mat3d deviatoric_stress_tensor = stress_tensor - (1.0 / stress_dimension_) * stress_tensor.trace() * Mat3d::Identity();
    Real stress_tensor_J2 = 0.5 * (deviatoric_stress_tensor.cwiseProduct(deviatoric_stress_tensor.transpose())).sum();
    if (-alpha_phi_ * stress_tensor_I1 + k_c_ < sqrt(stress_tensor_J2))
    {
        Real r = (-alpha_phi_ * stress_tensor_I1 + k_c_) / (sqrt(stress_tensor_J2) + TinyReal);
        stress_tensor = r * deviatoric_stress_tensor + (1.0 / stress_dimension_) * stress_tensor_I1 * Mat3d::Identity();
    }
    return stress_tensor;
}
//=================================================================================================//
J2Plasticity::J2Plasticity(Real rho0, Real c0, Real youngs_modulus, Real poisson_ratio,
                           Real yield_stress, Real hardening_modulus)
    : GeneralContinuum(rho0, c0, youngs_modulus, poisson_ratio),
      yield_stress_(yield_stress), hardening_modulus_(hardening_modulus),
      dv_hardening_factor_(nullptr)
{
    material_type_name_ = "J2Plasticity";
}
//=================================================================================================//
Matd J2Plasticity::ConstitutiveRelationShearStressWithHardening(Matd &velocity_gradient, Matd &shear_stress, Real &hardening_factor)
{
    Matd strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    Matd deviatoric_strain_rate = strain_rate - (1.0 / (Real)Dimensions) * strain_rate.trace() * Matd::Identity();
    Matd spin_rate = 0.5 * (velocity_gradient - velocity_gradient.transpose());
    Matd shear_stress_rate_elastic = 2.0 * G_ * deviatoric_strain_rate + shear_stress * (spin_rate.transpose()) + spin_rate * shear_stress;
    Real stress_tensor_J2 = 0.5 * (shear_stress.cwiseProduct(shear_stress.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor + yield_stress_);
    Real lambda_dot_ = 0;
    Matd g = Matd::Zero();
    if (f > TinyReal)
    {
        Real deviatoric_stress_times_strain_rate = (shear_stress.cwiseProduct(strain_rate)).sum();
        lambda_dot_ = deviatoric_stress_times_strain_rate / (sqrt(2.0 * stress_tensor_J2) * (1.0 + hardening_modulus_ / (3.0 * G_)));
        g = lambda_dot_ * (sqrt(2.0) * G_ * shear_stress / (sqrt(stress_tensor_J2)));
    }
    return shear_stress_rate_elastic - g;
}
//=================================================================================================//
Matd J2Plasticity::ReturnMappingShearStress(Matd &shear_stress, Real &hardening_factor)
{
    Real stress_tensor_J2 = 0.5 * (shear_stress.cwiseProduct(shear_stress.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor + yield_stress_);
    Real r = 1.0;
    if (f > TinyReal)
        r = (sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor + yield_stress_)) / (sqrt(2.0 * stress_tensor_J2) + TinyReal);
    return r * shear_stress;
}
//=================================================================================================//
Real J2Plasticity::ScalePenaltyForce(Matd &shear_stress, Real &hardening_factor)
{
    Real stress_tensor_J2 = 0.5 * (shear_stress.cwiseProduct(shear_stress.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor + yield_stress_);
    return (f > TinyReal) ? (sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor + yield_stress_)) / (sqrt(2.0 * stress_tensor_J2) + TinyReal) : 1.0;
}
//=================================================================================================//
Real J2Plasticity::HardeningFactorRate(const Matd &shear_stress, Real &hardening_factor)
{
    Real stress_tensor_J2 = 0.5 * (shear_stress.cwiseProduct(shear_stress.transpose())).sum();
    Real f = sqrt(2.0 * stress_tensor_J2) - sqrt_2_over_3_ * (hardening_modulus_ * hardening_factor + yield_stress_);
    return (f > TinyReal) ? 0.5 * f / (G_ + hardening_modulus_ / 3.0) : 0.0;
}
//=================================================================================================//
void J2Plasticity::initializeLocalParameters(BaseParticles *base_particles)
{
    GeneralContinuum::initializeLocalParameters(base_particles);
    dv_hardening_factor_ = base_particles->registerStateVariable<Real>("HardeningFactor");
    base_particles->addEvolvingVariable<Real>(dv_hardening_factor_);
}
//=================================================================================================//
} // namespace SPH