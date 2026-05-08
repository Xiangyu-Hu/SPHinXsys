#ifndef ELECTROMAGNETIC_APHI_GLOBAL_IMPLICIT_SOLVER_HPP
#define ELECTROMAGNETIC_APHI_GLOBAL_IMPLICIT_SOLVER_HPP

#include "electromagnetic_aphi_global_implicit_solver.h"
#include "electromagnetic_team7_aphi_frequency_dynamics.hpp"

namespace SPH
{
namespace electromagnetics
{
//=================================================================================================//
VectorPotentialFrequencyCoupledPreconditionerInner::
    VectorPotentialFrequencyCoupledPreconditionerInner(
        BaseInnerRelation &inner_relation,
        Real angular_frequency,
        const std::string &residual_real_name,
        const std::string &residual_imag_name,
        const std::string &search_real_name,
        const std::string &search_imag_name,
        Real sigma_relaxation_scaling,
        Real sigma_relaxation_floor,
        Real magnetic_diagonal_scaling,
        Real max_search_norm)
    : LocalDynamics(inner_relation.getSPHBody()),
      DataDelegateInner(inner_relation),
      angular_frequency_(angular_frequency),
      sigma_relaxation_scaling_(sigma_relaxation_scaling),
      sigma_relaxation_floor_(sigma_relaxation_floor),
      magnetic_diagonal_scaling_(magnetic_diagonal_scaling),
      max_search_norm_(max_search_norm),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      residual_real_(particles_->getVariableDataByName<Vecd>(residual_real_name)),
      residual_imag_(particles_->getVariableDataByName<Vecd>(residual_imag_name)),
      search_real_(particles_->getVariableDataByName<Vecd>(search_real_name)),
      search_imag_(particles_->getVariableDataByName<Vecd>(search_imag_name))
{
}
//=================================================================================================//
void VectorPotentialFrequencyCoupledPreconditionerInner::interaction(size_t index_i, Real dt)
{
    (void)dt;
    Real sigma_i = electrical_conductivity_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real conservative_diagonal =
        ComputeConservativeMagneticDiagonal(inner_neighborhood, index_i, Vol_,
                                            magnetic_reluctivity_, magnetic_diagonal_scaling_);
    Matd self_block =
        ComputeLocalMagneticSelfBlock(inner_neighborhood, index_i, Vol_,
                                      magnetic_reluctivity_, magnetic_diagonal_scaling_);
    Real fallback_operator_diagonal = conservative_diagonal + sigma_relaxation_floor_;
    Real omega_sigma = sigma_relaxation_scaling_ * sigma_i * angular_frequency_;

    Vecd delta_real = ZeroData<Vecd>::value;
    Vecd delta_imag = ZeroData<Vecd>::value;
    SolveCoupledMagneticSelfBlockOrFallback(
        self_block,
        residual_real_[index_i],
        residual_imag_[index_i],
        omega_sigma,
        fallback_operator_diagonal,
        delta_real,
        delta_imag);

    Real pair_norm = sqrt(delta_real.squaredNorm() + delta_imag.squaredNorm());
    if (!std::isfinite(pair_norm))
    {
        search_real_[index_i] = ZeroData<Vecd>::value;
        search_imag_[index_i] = ZeroData<Vecd>::value;
        return;
    }
    if (pair_norm > max_search_norm_)
    {
        Real scale = max_search_norm_ / (pair_norm + TinyReal);
        delta_real *= scale;
        delta_imag *= scale;
    }
    search_real_[index_i] = delta_real;
    search_imag_[index_i] = delta_imag;
}
//=================================================================================================//
FrequencyVectorPotentialLinearOperatorComplex::
    FrequencyVectorPotentialLinearOperatorComplex(
        SPHBody &sph_body,
        Real angular_frequency,
        const std::string &search_real_name,
        const std::string &search_imag_name,
        const std::string &curl_nu_b_real_name,
        const std::string &curl_nu_b_imag_name,
        const std::string &operator_action_real_name,
        const std::string &operator_action_imag_name)
    : LocalDynamics(sph_body),
      angular_frequency_(angular_frequency),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      search_real_(particles_->getVariableDataByName<Vecd>(search_real_name)),
      search_imag_(particles_->getVariableDataByName<Vecd>(search_imag_name)),
      curl_nu_b_real_(particles_->getVariableDataByName<Vecd>(curl_nu_b_real_name)),
      curl_nu_b_imag_(particles_->getVariableDataByName<Vecd>(curl_nu_b_imag_name)),
      operator_action_real_(particles_->getVariableDataByName<Vecd>(operator_action_real_name)),
      operator_action_imag_(particles_->getVariableDataByName<Vecd>(operator_action_imag_name))
{
}
//=================================================================================================//
void FrequencyVectorPotentialLinearOperatorComplex::update(size_t index_i, Real dt)
{
    (void)dt;
    Real sigma_i = electrical_conductivity_[index_i];
    Vecd action_real = curl_nu_b_real_[index_i] - sigma_i * angular_frequency_ * search_imag_[index_i];
    Vecd action_imag = curl_nu_b_imag_[index_i] + sigma_i * angular_frequency_ * search_real_[index_i];
    operator_action_real_[index_i] =
        std::isfinite(action_real.squaredNorm()) ? action_real : ZeroData<Vecd>::value;
    operator_action_imag_[index_i] =
        std::isfinite(action_imag.squaredNorm()) ? action_imag : ZeroData<Vecd>::value;
}
//=================================================================================================//

} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_GLOBAL_IMPLICIT_SOLVER_HPP
