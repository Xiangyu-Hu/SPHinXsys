/**
* @file 	diffusion_splitting_parameter.hpp
* @author	Bo Zhang and Xiangyu Hu
*/

#ifndef DIFFUSION_SPLITTING_PARAMETER_HPP
#define DIFFUSION_SPLITTING_PARAMETER_HPP

#include "diffusion_splitting_parameter.h"

namespace SPH
{	
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	ParameterSplittingByPDEInner<ParticlesType, VariableType>::
		ParameterSplittingByPDEInner(BaseInnerRelation& inner_relation, const std::string& variable_name):
		OptimizationBySplittingAlgorithmBase<ParticlesType, VariableType>(inner_relation, variable_name){};
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	ErrorAndParameters<VariableType> ParameterSplittingByPDEInner<ParticlesType, VariableType>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		VariableType& variable_i = this->variable_[index_i];
		ErrorAndParameters<VariableType> error_and_parameters;
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real& r_ij_ = inner_neighborhood.r_ij_[n];

			VariableType variable_derivative = variable_i + this->variable_[index_j];
			Real phi_ij = this->species_modified_[index_i] - this->species_recovery_[index_j];
			Real parameter_b = phi_ij * dW_ijV_j_ * dt / r_ij_;

			error_and_parameters.error_ -= variable_derivative * parameter_b;
			error_and_parameters.a_ += parameter_b;
			error_and_parameters.c_ += parameter_b * parameter_b;
		}
		error_and_parameters.a_ -= 1;
		error_and_parameters.error_ -= this->heat_source_[index_i] * dt;
		return error_and_parameters;
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	void ParameterSplittingByPDEInner<ParticlesType, VariableType>::
		 updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters)
	{
		Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
		VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
		this->parameter_recovery_[index_i] = this->variable_[index_i];
		this->variable_[index_i] += parameter_k * error_and_parameters.a_;
		if (this->variable_[index_i] < 0.1) { this->variable_[index_i] = 0.1; } //set lower bound

		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real& r_ij_ = inner_neighborhood.r_ij_[n];

			Real phi_ij = this->species_modified_[index_i] - this->species_recovery_[index_j];
			Real parameter_b = phi_ij * dW_ijV_j_ * dt / r_ij_;

			this->parameter_recovery_[index_j] = this->variable_[index_j];
			this->variable_[index_j] += parameter_k * parameter_b;
			if (this->variable_[index_j] < 0.1) { this->variable_[index_j] = 0.1; } //set lower bound
		}	
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	void ParameterSplittingByPDEInner<ParticlesType, VariableType>::interaction(size_t index_i, Real dt)
	{
		/* With Pseudo Time Step */
		this->splitting_index_[index_i] = 0;
		ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
		error_and_parameters.error_ = error_and_parameters.error_ - this->residual_T_local_[index_i];
		updateStatesByError(index_i, dt, error_and_parameters);
		this->residual_k_local_[index_i] = error_and_parameters.error_;

		ErrorAndParameters<VariableType> error_and_parameters_aft = computeErrorAndParameters(index_i, dt);
		error_and_parameters_aft.error_ = error_and_parameters_aft.error_ - this->residual_T_local_[index_i];
		this->residual_after_splitting_[index_i] = error_and_parameters_aft.error_;

		if (abs(this->residual_after_splitting_[index_i]) > abs(this->residual_k_local_[index_i]))
		{
			this->variable_[index_i] = this->parameter_recovery_[index_i];
			Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t& index_j = inner_neighborhood.j_[n];
				this->variable_[index_j] = this->parameter_recovery_[index_j];
			}

			ErrorAndParameters<VariableType> error_and_parameters_reverse = computeErrorAndParameters(index_i, -dt);
			error_and_parameters_reverse.error_ = error_and_parameters_reverse.error_ + this->residual_T_local_[index_i];
			updateStatesByError(index_i, -dt, error_and_parameters_reverse);
			this->residual_k_local_[index_i] = error_and_parameters_reverse.error_;

			ErrorAndParameters<VariableType> error_and_parameters_aft = computeErrorAndParameters(index_i, -dt);
			error_and_parameters_aft.error_ = error_and_parameters_aft.error_ + this->residual_T_local_[index_i];
			this->residual_after_splitting_[index_i] = error_and_parameters_aft.error_;

			if (abs(this->residual_after_splitting_[index_i]) > abs(this->residual_k_local_[index_i]))
			{
				this->splitting_index_[index_i] = 0;
				this->variable_[index_i] = this->parameter_recovery_[index_i];
				Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					size_t& index_j = inner_neighborhood.j_[n];
					this->variable_[index_j] = this->parameter_recovery_[index_j];
				}
			}
			else if (abs(this->residual_after_splitting_[index_i]) <= abs(this->residual_k_local_[index_i]))
			{
				this->splitting_index_[index_i] = 2;
			}
		}
		else if (abs(this->residual_after_splitting_[index_i]) <= abs(this->residual_k_local_[index_i]))
		{
			this->splitting_index_[index_i] = 1;
		}
	};
	//=================================================================================================//
	template <class ParticlesType, class ContactParticlesType, typename VariableType>
	ParameterSplittingByPDEWithBoundary<ParticlesType, ContactParticlesType, VariableType>::
		ParameterSplittingByPDEWithBoundary(ComplexRelation& complex_relation, const std::string& variable_name):
		ParameterSplittingByPDEInner<ParticlesType, VariableType>(complex_relation.getInnerRelation(), variable_name),
		DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>(complex_relation.getContactRelation())
	{
		boundary_heat_flux_.resize(this->contact_particles_.size());
		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			boundary_normal_vector_.push_back(&this->contact_particles_[k]->n_);
			boundary_species_.push_back(&(this->contact_particles_[k]->all_species_));
			boundary_heat_flux_[k] = this->contact_particles_[k]->template registerSharedVariable<Real>("HeatFlux");
		}
	}
	//=================================================================================================//
	template <class ParticlesType, class ContactParticlesType, typename VariableType>
	ErrorAndParameters<VariableType> ParameterSplittingByPDEWithBoundary<ParticlesType, ContactParticlesType, VariableType>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = ParameterSplittingByPDEInner<ParticlesType, VariableType>
			::computeErrorAndParameters(index_i, dt);

		VariableType& variable_i = this->variable_[index_i];
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& heat_flux_k = *(this->boundary_heat_flux_[k]);
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdVec<StdLargeVec<Real>>& species_k = *(boundary_species_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t& index_j = contact_neighborhood.j_[n];

				if (species_k[this->phi_][index_j] > 0.0)
				{
					VariableType variable_derivative = variable_i;
					Real phi_ij = this->species_modified_[index_i] - species_k[this->phi_][index_j];
					Real parameter_b = 2.0 * phi_ij * contact_neighborhood.dW_ijV_j_[n] * dt / contact_neighborhood.r_ij_[n];

					error_and_parameters.error_ -= variable_derivative * parameter_b;
					error_and_parameters.a_ += parameter_b;
				}

				if (heat_flux_k[index_j] != 0.0)
				{
					Vecd n_ij = this->normal_vector_[index_i] - normal_vector_k[index_j];
					error_and_parameters.error_ -= heat_flux_k[index_j] * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n].dot(n_ij) * dt;
				}
			}
		}
		return error_and_parameters;
	}
	//=================================================================================================//
	template <typename ParameterSplittingType, typename BaseBodyRelationType, typename VariableType>
	UpdateParameterPDEResidual<ParameterSplittingType, BaseBodyRelationType, VariableType>::
		UpdateParameterPDEResidual(BaseBodyRelationType& body_relation, const std::string& variable_name):
		ParameterSplittingType(body_relation, variable_name) {};
	//=================================================================================================//
	template <typename ParameterSplittingType, typename BaseBodyRelationType, typename VariableType>
	void UpdateParameterPDEResidual<ParameterSplittingType, BaseBodyRelationType, VariableType>::
		interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = this->computeErrorAndParameters(index_i, dt);
		error_and_parameters.error_ = error_and_parameters.error_ - this->residual_T_local_[index_i];
		this->residual_k_global_[index_i] = error_and_parameters.error_;
	}
	//=================================================================================================//
}
#endif //DIFFUSION_SPLITTING_PARAMETER_HPP
