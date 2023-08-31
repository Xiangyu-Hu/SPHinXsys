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
	template <class ParticleType, typename VariableType>
	OptimizationBySplittingAlgorithmBase<ParticleType, VariableType>::
		OptimizationBySplittingAlgorithmBase(BaseInnerRelation& inner_relation, const std::string &variable_name) :
		BaseDiffusionRelaxation<ParticleType>(inner_relation.getSPHBobdy()),
		DataDelegateInner<ParticleType, DataDelegateEmptyBase>(inner_relation),
		Vol_(this->particles_->Vol_),
		mass_(this->particles_->mass_),
		heat_flux_(this->particles_->heat_flux_),
		heat_source_(this->particles_->heat_source_),
		variable_(*this->particles_->template getVariableByName<VariableType>(variable_name))
	{
		*this->particles_->registerSharedVariable(splitting_index_, "splitting_index");
		*this->particles_->addVariableToWrite<Real>("splitting_index");

		*this->particles_->registerSharedVariable(species_modified_, "species_modified");
		*this->partivles_->addVariableToWrite<Real>("species_modified");

		*this->particles_->registerSharedVariable(species_recovery_, "species_recovery");
		*this->particles_->addVariableToWrite<Real>("species_recovery");

		*this->particles_->registerSharedVariable(parameter_recovery_, "parameter_recovery");
		*this->particles_->addVariableToWrite<Real>("parameter_recovery");

		*this->particles_->registerSharedVariable(eta_regularization_, "eta_regularization");
		*this->particles_->addVariableToWrite<Real>("eta_regularization");

		*this->particles_->registerSharedVariable(normal_distance_, "normal_distance");
		*this->particles_->addVariableToWrite<Real>("normal_distance");

		*this->particles_->registerSharedVariable(normal_vector_, "normal_vector");
		*this->particles_->addVariableToWrite<Real>("normal_vector");

		*this->particles_->registerSharedVariable(residual_T_local_, "residual_T_local");
		*this->particles_->addVariableToWrite<Real>("residual_T_local");

		*this->particles_->registerSharedVariable(residual_T_global_, "residual_T_global");
		*this->particles_->addVariableToWrite<Real>("residual_T_global");

		*this->particles_->registerSharedVariable(residual_k_local_, "residual_k_local");
		*this->particles_->addVariableToWrite<Real>("residual_k_local");

		*this->particles_->registerSharedVariable(residual_k_global_, "residual_k_global");
		*this->particles_->addVariableToWrite<Real>("residual_k_global");

		*this->particles_->registerSharedVariable(variation_local_, "variation_local");
		*this->particles_->addVariableToWrite<Real>("variation_local");

		*this->particles_->registerSharedVariable(variation_global_, "variation_global");
		*this->particles_->addVariableToWrite<Real>("variation_global");

		phi_ = this->particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class ParticleType, typename VariableType>
	RegularizationByDiffusionAnalogy<ParticleType, VariableType>::
		RegularizationByDiffusionAnalogy(BaseInnerRelation& inner_relation, const std::string& variable_name, 
			                             Real initial_eta, Real variation) :
		OptimizationBySplittingAlgorithmBase<ParticleType, VariableType>(inner_relation, variable_name),
		initial_eta_(initial_eta), maximum_variation_(variation), averaged_variation_(variation) {}
	//=================================================================================================//
	template <class ParticleType, typename VariableType>
	ErrorAndParameters<VariableType> RegularizationByDiffusionAnalogy<ParticleType, VariableType>::
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
	template <class ParticleType, typename VariableType>
	void RegularizationByDiffusionAnalogy<ParticleType, VariableType>::
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
	template <class ParticleType, typename VariableType>
	void RegularizationByDiffusionAnalogy<ParticleType, VariableType>::interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = computeVariationAndParameters(index_i, dt);
		updateStatesByVariation(index_i, dt, error_and_parameters);
		this->variation_local_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
	}
	//=================================================================================================//
	template <class ParticleType, typename VariableType>
	UpdateRegularizationVariation<ParticleType, VariableType>::
		UpdateRegularizationVariation(BaseInnerRelation& inner_relation, const std::string& variable_name) :
		OptimizationBySplittingAlgorithmBase<ParticleType, VariableType>(inner_relation, variable_name) {};
	//=================================================================================================//
	template <class ParticleType, typename VariableType>
	ErrorAndParameters<VariableType> UpdateRegularizationVariation<ParticleType, VariableType>::
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
	template <class ParticleType, typename VariableType>
	void UpdateRegularizationVariation<ParticleType, VariableType>::interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = this->computeVariationAndParameters(index_i, dt);
		this->variation_global_[index_i] = error_and_parameters.error_ / dt / this->eta_regularization_[index_i];
	}
	//=================================================================================================//
};

#endif DIFFUSION_SPLITTING_BASE_H