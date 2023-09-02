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
	template <class ParticlesType, typename VariableType>
	OptimizationBySplittingAlgorithmBase<ParticlesType, VariableType>::
		OptimizationBySplittingAlgorithmBase(BaseInnerRelation& inner_relation, const std::string &variable_name) :
		BaseDiffusionRelaxation<ParticleType>(inner_relation.getSPHBobdy()),
		DataDelegateInner<ParticlesType, DataDelegateEmptyBase>(inner_relation),
		Vol_(this->particles_->Vol_), mass_(this->particles_->mass_), normal_vector_(this->particles_->n_),
		heat_flux_(this->particles_->heat_flux_), heat_source_(this->particles_->heat_source_),
		variable_(*this->particles_->template getVariableByName<VariableType>(variable_name))
	{
        *this->particles_->registerSharedVariable(heat_flux_, "HeatFlux");
        *this->particles_->addVariableToWrite<Real>("HeatFlux");

		*this->particles_->registerSharedVariable(heat_flux_, "HeatSource");
        *this->particles_->addVariableToWrite<Real>("HeatSource");

		*this->particles_->registerSharedVariable(splitting_index_, "SplittingIndex");
		*this->particles_->addVariableToWrite<Real>("SplittingIndex");

		*this->particles_->registerSharedVariable(species_modified_, "SpeciesModified");
		*this->partivles_->addVariableToWrite<Real>("SpeciesModified");

		*this->particles_->registerSharedVariable(species_recovery_, "SpeciesRecovery");
		*this->particles_->addVariableToWrite<Real>("SpeciesRecovery");

		*this->particles_->registerSharedVariable(parameter_recovery_, "ParameterRecovery");
		*this->particles_->addVariableToWrite<Real>("ParameterRecovery");

		*this->particles_->registerSharedVariable(eta_regularization_, "EtaRegularization");
		*this->particles_->addVariableToWrite<Real>("EtaRegularization");

		*this->particles_->registerSharedVariable(residual_T_local_, "ResidualTLocal");
		*this->particles_->addVariableToWrite<Real>("ResidualTLocal");

		*this->particles_->registerSharedVariable(residual_T_global_, "ResidualTGlobal");
		*this->particles_->addVariableToWrite<Real>("ResidualTGlobal");

		*this->particles_->registerSharedVariable(residual_k_local_, "ResidualKLocal");
		*this->particles_->addVariableToWrite<Real>("ResidualKLocal");

		*this->particles_->registerSharedVariable(residual_k_global_, "ResidualKGlobal");
		*this->particles_->addVariableToWrite<Real>("ResidualKGlobal");

		*this->particles_->registerSharedVariable(variation_local_, "VariationLocal");
		*this->particles_->addVariableToWrite<Real>("VariationLocal");

		*this->particles_->registerSharedVariable(variation_global_, "VariationGlobal");
		*this->particles_->addVariableToWrite<Real>("VariationGlobal");

		phi_ = this->particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	RegularizationByDiffusionAnalogy<ParticlesType, VariableType>::
		RegularizationByDiffusionAnalogy(BaseInnerRelation& inner_relation, const std::string& variable_name, 
			                             Real initial_eta, Real variation) :
		OptimizationBySplittingAlgorithmBase<ParticlesType, VariableType>(inner_relation, variable_name),
		initial_eta_(initial_eta), maximum_variation_(variation), averaged_variation_(variation) {}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	ErrorAndParameters<VariableType> RegularizationByDiffusionAnalogy<ParticlesType, VariableType>::
		computeVariationAndParameters(size_t index_i, Real dt)
	{
		Real Vol_i = this->Vol_[index_i];
		Real mass_i = this->mass_[index_i];
		VariableType& variable_i = this->variable_[index_i];
		ErrorAndParameters<VariableType> error_and_parameters;
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real& r_ij_ = inner_neighborhood.r_ij_[n];

			//this->eta_regularization_[index_i] = initial_eta_ * abs(this->variation_local_[index_i] + TinyReal) / averaged_variation_;
			//this->eta_regularization_[index_i] = initial_eta_ * abs(this->variation_local_[index_i] + TinyReal) / abs(maximum_variation_);
			this->eta_regularization_[index_i] = initial_eta_; //uniform coefficient.

			VariableType variable_derivative = (variable_i - this->variable_[index_j]);
			Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j_ * Vol_i * dt / r_ij_;

			error_and_parameters.error_ -= variable_derivative * parameter_b;
			error_and_parameters.a_ += parameter_b;
			error_and_parameters.c_ += parameter_b * parameter_b;
		}
		error_and_parameters.a_ -= mass_i;
		return error_and_parameters;
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	void RegularizationByDiffusionAnalogy<ParticlesType, VariableType>::
		updateStatesByVariation(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters)
	{
		Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
		VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
		this->variable_[index_i] += parameter_k * error_and_parameters.a_;
		if (this->variable_[index_i] < 0.1) { this->variable_[index_i] = 0.01; } //set lower bound.

		Real Vol_i = this->Vol_[index_i];
		VariableType& variable_i = this->variable_[index_i];
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real& r_ij_ = inner_neighborhood.r_ij_[n];

			Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j_ * Vol_i * dt / r_ij_;

			//predicted quantity at particle j
			VariableType variable_j = this->variable_[index_j] - parameter_k * parameter_b;
			VariableType variable_derivative = (variable_i - variable_j);

			//exchange in conservation form
			this->variable_[index_j] -= variable_derivative * parameter_b / this->mass_[index_j];
			if (this->variable_[index_j] < 0.1) { this->variable_[index_j] = 0.01; } //set lower bound.
		}
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	void RegularizationByDiffusionAnalogy<ParticlesType, VariableType>::interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = computeVariationAndParameters(index_i, dt);
		updateStatesByVariation(index_i, dt, error_and_parameters);
		this->variation_local_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	UpdateRegularizationVariation<ParticlesType, VariableType>::
		UpdateRegularizationVariation(BaseInnerRelation& inner_relation, const std::string& variable_name) :
		OptimizationBySplittingAlgorithmBase<ParticleType, VariableType>(inner_relation, variable_name) {};
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	ErrorAndParameters<VariableType> UpdateRegularizationVariation<ParticlesType, VariableType>::
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
			Real &dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real &r_ij_ = inner_neighborhood.r_ij_[n];

			VariableType variable_derivative = (variable_i - this->variable_[index_j]);
			Real parameter_b = 2.0 * this->eta_regularization_[index_i] * dW_ijV_j_ * Vol_i * dt / r_ij_;
		                  	   
			error_and_parameters.error_ -= variable_derivative * parameter_b;
			error_and_parameters.a_ += parameter_b;
			error_and_parameters.c_ += parameter_b * parameter_b;
		}
		error_and_parameters.a_ -= mass_i;
		return error_and_parameters;
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	void UpdateRegularizationVariation<ParticlesType, VariableType>::interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = this->computeVariationAndParameters(index_i, dt);
		this->variation_global_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
	}
	//=================================================================================================//
};

#endif DIFFUSION_SPLITTING_BASE_H