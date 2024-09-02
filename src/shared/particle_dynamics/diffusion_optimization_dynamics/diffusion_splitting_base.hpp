/**
 * @file 	diffusion_splitting_base.hpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_SPLITTING_BASE_HPP
#define DIFFUSION_SPLITTING_BASE_HPP

#include "diffusion_splitting_base.h"

namespace SPH
{
//=================================================================================================//
template <typename VariableType>
OptimizationBySplittingAlgorithmBase<VariableType>::
    OptimizationBySplittingAlgorithmBase(BaseInnerRelation &inner_relation, const std::string &variable_name)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      diffusion_(DynamicCast<LocalIsotropicDiffusion>(this, sph_body_.getBaseMaterial())),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(*this->particles_->template getVariableDataByName<Real>("Mass")),
      normal_vector_(*this->particles_->template getVariableDataByName<Vecd>("NormalDirection")),
      variable_(*(this->particles_->template registerSharedVariable<VariableType>(variable_name))),
      heat_flux_(*(this->particles_->template registerSharedVariable<Real>("HeatFlux"))),
      heat_source_(*(this->particles_->template registerSharedVariable<Real>("HeatSource"))),
      splitting_index_(*(this->particles_->template registerSharedVariable<int>("SplittingIndex"))),
      species_modified_(*(this->particles_->template registerSharedVariable<Real>("SpeciesModified"))),
      species_recovery_(*(this->particles_->template registerSharedVariable<Real>("SpeciesRecovery"))),
      parameter_recovery_(*(this->particles_->template registerSharedVariable<Real>("ParameterRecovery"))),
      eta_regularization_(*(this->particles_->template registerSharedVariable<Real>("EtaRegularization"))),
      residual_T_local_(*(this->particles_->template registerSharedVariable<Real>("ResidualTLocal"))),
      residual_T_global_(*(this->particles_->template registerSharedVariable<Real>("ResidualTGlobal"))),
      residual_k_local_(*(this->particles_->template registerSharedVariable<Real>("ResidualKLocal"))),
      residual_k_global_(*(this->particles_->template registerSharedVariable<Real>("ResidualKGlobal"))),
      variation_local_(*(this->particles_->template registerSharedVariable<Real>("VariationLocal"))),
      variation_global_(*(this->particles_->template registerSharedVariable<Real>("VariationGlobal"))),
      residual_after_splitting_(*(this->particles_->template registerSharedVariable<Real>("ResidualAfterSplitting")))
{
    this->particles_->template addVariableToWrite<Real>("HeatFlux");
    this->particles_->template addVariableToWrite<Real>("HeatSource");
    this->particles_->template addVariableToWrite<int>("SplittingIndex");
    this->particles_->template addVariableToWrite<Real>("SpeciesModified");
    this->particles_->template addVariableToWrite<Real>("SpeciesRecovery");
    this->particles_->template addVariableToWrite<Real>("ParameterRecovery");
    this->particles_->template addVariableToWrite<Real>("EtaRegularization");
    this->particles_->template addVariableToWrite<Real>("ResidualTLocal");
    this->particles_->template addVariableToWrite<Real>("ResidualTGlobal");
    this->particles_->template addVariableToWrite<Real>("ResidualKLocal");
    this->particles_->template addVariableToWrite<Real>("ResidualKGlobal");
    this->particles_->template addVariableToWrite<Real>("VariationLocal");
    this->particles_->template addVariableToWrite<Real>("VariationGlobal");
    this->particles_->template addVariableToWrite<Real>("ResidualAfterSplitting");
    this->particles_->template addVariableToRestart<Real>("ThermalConductivity");
}
//=================================================================================================//
template <typename VariableType>
RegularizationByDiffusionAnalogy<VariableType>::
    RegularizationByDiffusionAnalogy(BaseInnerRelation &inner_relation, const std::string &variable_name,
                                     Real initial_eta, Real variation)
    : OptimizationBySplittingAlgorithmBase<VariableType>(inner_relation, variable_name),
      initial_eta_(initial_eta), maximum_variation_(variation), averaged_variation_(variation) {}
//=================================================================================================//
template <typename VariableType>
ErrorAndParameters<VariableType> RegularizationByDiffusionAnalogy<VariableType>::
    computeVariationAndParameters(size_t index_i, Real dt)
{
    Real Vol_i = this->Vol_[index_i];
    Real mass_i = this->mass_[index_i];
    VariableType &variable_i = this->variable_[index_i];
    ErrorAndParameters<VariableType> error_and_parameters;
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij_ = inner_neighborhood.r_ij_[n];

        // this->eta_regularization_[index_i] = initial_eta_ * abs(this->variation_local_[index_i] + TinyReal) / averaged_variation_;
        // this->eta_regularization_[index_i] = initial_eta_ * abs(this->variation_local_[index_i] + TinyReal) / abs(maximum_variation_);
        this->eta_regularization_[index_i] = initial_eta_; // uniform coefficient.

        VariableType variable_derivative = (variable_i - this->variable_[index_j]);
        Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j * Vol_i * dt / r_ij_;

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= mass_i;
    return error_and_parameters;
}
//=================================================================================================//
template <typename VariableType>
void RegularizationByDiffusionAnalogy<VariableType>::
    updateStatesByVariation(size_t index_i, Real dt, const ErrorAndParameters<VariableType> &error_and_parameters)
{
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    this->variable_[index_i] += parameter_k * error_and_parameters.a_;
    if (this->variable_[index_i] < 0.1)
    {
        this->variable_[index_i] = 0.1;
    } // set lower bound.

    Real Vol_i = this->Vol_[index_i];
    VariableType &variable_i = this->variable_[index_i];
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];

        Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j * Vol_i * dt / r_ij;

        // predicted quantity at particle j
        VariableType variable_j = this->variable_[index_j] - parameter_k * parameter_b;
        VariableType variable_derivative = (variable_i - variable_j);

        // exchange in conservation form
        this->variable_[index_j] -= variable_derivative * parameter_b / this->mass_[index_j];
        if (this->variable_[index_j] < 0.1)
        {
            this->variable_[index_j] = 0.1;
        } // set lower bound.
    }
}
//=================================================================================================//
template <typename VariableType>
void RegularizationByDiffusionAnalogy<VariableType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters = computeVariationAndParameters(index_i, dt);
    updateStatesByVariation(index_i, dt, error_and_parameters);
    this->variation_local_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
}
//=================================================================================================//
template <typename VariableType>
UpdateRegularizationVariation<VariableType>::
    UpdateRegularizationVariation(BaseInnerRelation &inner_relation, const std::string &variable_name)
    : OptimizationBySplittingAlgorithmBase<VariableType>(inner_relation, variable_name){};
//=================================================================================================//
template <typename VariableType>
ErrorAndParameters<VariableType> UpdateRegularizationVariation<VariableType>::
    computeVariationAndParameters(size_t index_i, Real dt)
{
    Real Vol_i = this->Vol_[index_i];
    Real mass_i = this->mass_[index_i];
    VariableType &variable_i = this->variable_[index_i];
    ErrorAndParameters<VariableType> error_and_parameters;

    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];

        VariableType variable_derivative = (variable_i - this->variable_[index_j]);
        Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j * Vol_i * dt / r_ij;

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= mass_i;
    return error_and_parameters;
}
//=================================================================================================//
template <typename VariableType>
void UpdateRegularizationVariation<VariableType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<VariableType> error_and_parameters = this->computeVariationAndParameters(index_i, dt);
    this->variation_global_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
}
//=================================================================================================//
}; // namespace SPH

#endif // DIFFUSION_SPLITTING_BASE_H