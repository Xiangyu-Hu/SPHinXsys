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
template <typename DataType>
OptimizationBySplittingAlgorithmBase<DataType>::
    OptimizationBySplittingAlgorithmBase(BaseInnerRelation &inner_relation, const std::string &variable_name)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      diffusion_(DynamicCast<LocalIsotropicDiffusion>(this, sph_body_.getBaseMaterial())),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      normal_vector_(this->particles_->template getVariableDataByName<Vecd>("NormalDirection")),
      variable_(this->particles_->template registerStateVariable<DataType>(variable_name)),
      heat_flux_(this->particles_->template registerStateVariable<Real>("HeatFlux")),
      heat_source_(this->particles_->template registerStateVariable<Real>("HeatSource")),
      splitting_index_(this->particles_->template registerStateVariable<int>("SplittingIndex")),
      species_modified_(this->particles_->template registerStateVariable<Real>("SpeciesModified")),
      species_recovery_(this->particles_->template registerStateVariable<Real>("SpeciesRecovery")),
      parameter_recovery_(this->particles_->template registerStateVariable<Real>("ParameterRecovery")),
      eta_regularization_(this->particles_->template registerStateVariable<Real>("EtaRegularization")),
      residual_T_local_(this->particles_->template registerStateVariable<Real>("ResidualTLocal")),
      residual_T_global_(this->particles_->template registerStateVariable<Real>("ResidualTGlobal")),
      residual_k_local_(this->particles_->template registerStateVariable<Real>("ResidualKLocal")),
      residual_k_global_(this->particles_->template registerStateVariable<Real>("ResidualKGlobal")),
      variation_local_(this->particles_->template registerStateVariable<Real>("VariationLocal")),
      variation_global_(this->particles_->template registerStateVariable<Real>("VariationGlobal")),
      residual_after_splitting_(this->particles_->template registerStateVariable<Real>("ResidualAfterSplitting"))
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
    this->particles_->template addEvolvingVariable<Real>("ThermalConductivity");
}
//=================================================================================================//
template <typename DataType>
RegularizationByDiffusionAnalogy<DataType>::
    RegularizationByDiffusionAnalogy(BaseInnerRelation &inner_relation, const std::string &variable_name,
                                     Real initial_eta, Real variation)
    : OptimizationBySplittingAlgorithmBase<DataType>(inner_relation, variable_name),
      initial_eta_(initial_eta), maximum_variation_(variation), averaged_variation_(variation) {}
//=================================================================================================//
template <typename DataType>
ErrorAndParameters<DataType> RegularizationByDiffusionAnalogy<DataType>::
    computeVariationAndParameters(size_t index_i, Real dt)
{
    Real Vol_i = this->Vol_[index_i];
    Real mass_i = this->mass_[index_i];
    DataType &variable_i = this->variable_[index_i];
    ErrorAndParameters<DataType> error_and_parameters;
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij_ = inner_neighborhood.r_ij_[n];

        // this->eta_regularization_[index_i] = initial_eta_ * abs(this->variation_local_[index_i] + TinyReal) / averaged_variation_;
        // this->eta_regularization_[index_i] = initial_eta_ * abs(this->variation_local_[index_i] + TinyReal) / abs(maximum_variation_);
        this->eta_regularization_[index_i] = initial_eta_; // uniform coefficient.

        DataType variable_derivative = (variable_i - this->variable_[index_j]);
        Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j * Vol_i * dt / r_ij_;

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= mass_i;
    return error_and_parameters;
}
//=================================================================================================//
template <typename DataType>
void RegularizationByDiffusionAnalogy<DataType>::
    updateStatesByVariation(size_t index_i, Real dt, const ErrorAndParameters<DataType> &error_and_parameters)
{
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    DataType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    this->variable_[index_i] += parameter_k * error_and_parameters.a_;
    if (this->variable_[index_i] < 0.1)
    {
        this->variable_[index_i] = 0.1;
    } // set lower bound.

    Real Vol_i = this->Vol_[index_i];
    DataType &variable_i = this->variable_[index_i];
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];

        Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j * Vol_i * dt / r_ij;

        // predicted quantity at particle j
        DataType variable_j = this->variable_[index_j] - parameter_k * parameter_b;
        DataType variable_derivative = (variable_i - variable_j);

        // exchange in conservation form
        this->variable_[index_j] -= variable_derivative * parameter_b / this->mass_[index_j];
        if (this->variable_[index_j] < 0.1)
        {
            this->variable_[index_j] = 0.1;
        } // set lower bound.
    }
}
//=================================================================================================//
template <typename DataType>
void RegularizationByDiffusionAnalogy<DataType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<DataType> error_and_parameters = computeVariationAndParameters(index_i, dt);
    updateStatesByVariation(index_i, dt, error_and_parameters);
    this->variation_local_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
}
//=================================================================================================//
template <typename DataType>
UpdateRegularizationVariation<DataType>::
    UpdateRegularizationVariation(BaseInnerRelation &inner_relation, const std::string &variable_name)
    : OptimizationBySplittingAlgorithmBase<DataType>(inner_relation, variable_name){};
//=================================================================================================//
template <typename DataType>
ErrorAndParameters<DataType> UpdateRegularizationVariation<DataType>::
    computeVariationAndParameters(size_t index_i, Real dt)
{
    Real Vol_i = this->Vol_[index_i];
    Real mass_i = this->mass_[index_i];
    DataType &variable_i = this->variable_[index_i];
    ErrorAndParameters<DataType> error_and_parameters;

    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * this->Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];

        DataType variable_derivative = (variable_i - this->variable_[index_j]);
        Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j * Vol_i * dt / r_ij;

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= mass_i;
    return error_and_parameters;
}
//=================================================================================================//
template <typename DataType>
void UpdateRegularizationVariation<DataType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<DataType> error_and_parameters = this->computeVariationAndParameters(index_i, dt);
    this->variation_global_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
}
//=================================================================================================//
}; // namespace SPH

#endif // DIFFUSION_SPLITTING_BASE_H