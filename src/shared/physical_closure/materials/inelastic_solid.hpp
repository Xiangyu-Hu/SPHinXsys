#ifndef INELASTIC_SOLID_HPP
#define INELASTIC_SOLID_HPP

#include "inelastic_solid.h"

namespace SPH
{
//=================================================================================================//
template <typename ExecutionPolicy>
HardeningPlasticSolid::ConstituteKernel::ConstituteKernel(
    const ExecutionPolicy &ex_policy, HardeningPlasticSolid &encloser)
    : rho0_(encloser.ReferenceDensity()), K0_(encloser.BulkModulus()), G0_(encloser.ShearModulus()),
      c0_(encloser.ReferenceSoundSpeed()), cs0_(encloser.ShearWaveSpeed()),
      hardening_modulus_(encloser.HardeningModulus()), yield_stress_(encloser.YieldStress()),
      inverse_plastic_strain_(encloser.dv_inverse_plastic_strain_->DelegatedData(ex_policy)),
      hardening_parameter_(encloser.dv_hardening_parameter_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Matd HardeningPlasticSolid::ConstituteKernel::ElasticLeftCauchy(
    const Matd &F, size_t index_i, Real dt)
{
    Matd be = F * inverse_plastic_strain_[index_i] * F.transpose();
    Matd normalized_be = be * pow(be.determinant(), -OneOverDimensions);
    Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
    Matd deviatoric_Kirchhoff = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
    Real deviatoric_Kirchhoff_norm = deviatoric_Kirchhoff.norm();
    Real trial_function = deviatoric_Kirchhoff_norm -
                          sqrt_2_over_3_ * (hardening_modulus_ * hardening_parameter_[index_i] + yield_stress_);
    if (trial_function > 0.0)
    {
        Real renormalized_shear_modulus = normalized_be_isentropic * G0_;
        Real relax_increment = 0.5 * trial_function / (renormalized_shear_modulus + hardening_modulus_ / 3.0);
        hardening_parameter_[index_i] += sqrt_2_over_3_ * relax_increment;
        deviatoric_Kirchhoff -= 2.0 * renormalized_shear_modulus * relax_increment * deviatoric_Kirchhoff / deviatoric_Kirchhoff_norm;
        Matd relaxed_be = deviatoric_Kirchhoff / G0_ + normalized_be_isentropic * Matd::Identity();
        normalized_be = relaxed_be * pow(relaxed_be.determinant(), -OneOverDimensions);
    }
    Matd inverse_F = F.inverse();
    Matd inverse_F_T = inverse_F.transpose();
    inverse_plastic_strain_[index_i] = inverse_F * normalized_be * inverse_F_T;

    return normalized_be;
}
//=================================================================================================//
inline Real HardeningPlasticSolid::ConstituteKernel::VolumetricKirchhoff(Real J)
{
    return 0.5 * K0_ * (J * J - 1);
}
//=================================================================================================//
inline Matd HardeningPlasticSolid::ConstituteKernel::DeviatoricKirchhoff(const Matd &deviatoric_be)
{
    return G0_ * deviatoric_be;
}
//=================================================================================================//
template <typename ScalingType>
Matd HardeningPlasticSolid::ConstituteKernel::NumericalDampingLeftCauchy(
    const Matd &deformation, const Matd &deformation_rate, const ScalingType &scaling, size_t index_i)
{
    Matd strain_rate = 0.5 * (deformation_rate * deformation.transpose() +
                              deformation * deformation_rate.transpose());
    Matd normal_rate = getDiagonal(strain_rate);
    return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * scaling;
}
//=================================================================================================//
} // namespace SPH
#endif // INELASTIC_SOLID_HPP