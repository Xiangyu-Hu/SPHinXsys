#include "johnson_cook_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
JohnsonCookIntegration1stHalf::JohnsonCookIntegration1stHalf(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation), johnson_cook_solid_(DynamicCast<JohnsonCookSolid>(this, elastic_solid_)),
      material_active_(particles_->getVariableDataByName<int>("MaterialActive")) {}
//=================================================================================================//
void JohnsonCookIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * Real(0.5);
    if (!material_active_[index_i])
    {
        stress_PK1_B_[index_i].setZero();
        return;
    }
    F_[index_i] += dF_dt_[index_i] * dt * Real(0.5);
    const Matd sigma = johnson_cook_solid_.UpdateStress(F_[index_i], index_i, dt);
    const Real J = F_[index_i].determinant();
    rho_[index_i] = rho0_ / J;
    stress_PK1_B_[index_i] = J * sigma * F_[index_i].inverse().transpose() * B_[index_i].transpose();
}
//=================================================================================================//
void JohnsonCookIntegration1stHalf::interaction(size_t index_i, Real)
{
    Vecd force = Vecd::Zero();
    if (material_active_[index_i])
    {
        const Neighborhood &neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != neighborhood.current_size_; ++n)
        {
            const size_t index_j = neighborhood.j_[n];
            if (!material_active_[index_j])
                continue;
            const Real inverse_distance = Dimensions / neighborhood.r_ij_[n];
            const Real strain_rate = inverse_distance * inverse_distance *
                                     (pos_[index_i] - pos_[index_j]).dot(vel_[index_i] - vel_[index_j]);
            const Real weight = neighborhood.W_ij_[n] * inv_W0_;
            const Matd numerical_stress = Real(0.5) * (F_[index_i] + F_[index_j]) *
                                          elastic_solid_.PairNumericalDamping(strain_rate, smoothing_length_);
            force += mass_[index_i] * inv_rho0_ * neighborhood.dW_ij_[n] * Vol_[index_j] *
                     (stress_PK1_B_[index_i] + stress_PK1_B_[index_j] + numerical_dissipation_factor_ * weight * numerical_stress) * neighborhood.e_ij_[n];
        }
    }
    force_[index_i] = force;
}
//=================================================================================================//
JohnsonCookIntegration2ndHalf::JohnsonCookIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : Integration2ndHalf(inner_relation), material_active_(particles_->getVariableDataByName<int>("MaterialActive")) {}
//=================================================================================================//
void JohnsonCookIntegration2ndHalf::interaction(size_t index_i, Real)
{
    Matd deformation_rate = Matd::Zero();
    if (material_active_[index_i])
    {
        const Neighborhood &neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != neighborhood.current_size_; ++n)
        {
            const size_t index_j = neighborhood.j_[n];
            if (!material_active_[index_j])
                continue;
            const Vecd gradient = neighborhood.dW_ij_[n] * Vol_[index_j] * neighborhood.e_ij_[n];
            deformation_rate -= (vel_[index_i] - vel_[index_j]) * gradient.transpose();
        }
    }
    dF_dt_[index_i] = deformation_rate * B_[index_i];
}
} // namespace solid_dynamics
} // namespace SPH
