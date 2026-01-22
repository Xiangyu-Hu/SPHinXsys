#ifndef SHEAR_INTEGRATION_HPP
#define SHEAR_INTEGRATION_HPP

#include "shear_integration.h"

#include "base_particles.hpp"

namespace SPH
{
namespace continuum_dynamics
{
//====================================================================================//
template <class MaterialType, typename... Parameters>
void ShearForce<Inner<WithInitialization, MaterialType, Parameters...>>::
    InitializeKernel::update(size_t index_i, Real dt)
{
    Matd shear_stress_rate = J2_plasticity_.ConstitutiveRelationShearStressWithHardening(
        velocity_gradient_[index_i], shear_stress_[index_i], hardening_factor_[index_i]);
    Matd shear_stress_try = shear_stress_[index_i] + shear_stress_rate * dt;
    Real hardening_factor_increment = J2_plasticity_.HardeningFactorRate(shear_stress_try, hardening_factor_[index_i]);
    hardening_factor_[index_i] += sqrt(2.0 / 3.0) * hardening_factor_increment;
    scale_penalty_force_[index_i] = xi_ * J2_plasticity_.ScalePenaltyForce(shear_stress_try, hardening_factor_[index_i]);
    shear_stress_[index_i] = J2_plasticity_.ReturnMappingShearStress(shear_stress_try, hardening_factor_[index_i]);
    Matd strain_rate = 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_i].transpose());
    strain_tensor_[index_i] += strain_rate * dt;
}
//====================================================================================//
void ShearStressRelaxationHourglassControl1stHalf::interaction(size_t index_i, Real dt)
{
    Matd velocity_gradient = Matd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Vecd vel_i = vel_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Matd velocity_gradient_ij = -(vel_i - vel_[index_j]) * (B_[index_i] * e_ij * dW_ijV_j).transpose();
        velocity_gradient += velocity_gradient_ij;
    }
    velocity_gradient_[index_i] = velocity_gradient;
}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
#endif // SHEAR_INTEGRATION_HPP