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
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	ParameterSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		ParameterSplittingByPDEInner(BaseInnerRelation& inner_relation, const std::string& variable_name):
		OptimizationBySplittingAlgorithmBase<BaseParticlesType, BaseMaterialType, 
		                                     VariableType, NUM_SPECIES>(inner_relation, variable_name) {};
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> ParameterSplittingByPDEInner<BaseParticlesType, BaseMaterialType, 
		VariableType, NUM_SPECIES>::computeErrorAndParameters(size_t index_i, Real dt)
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
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ParameterSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
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
			this->variable_[index_j] += parameter_k * parameter_b; //the particle j is also regularized
			if (this->variable_[index_j] < 0.1) { this->variable_[index_j] = 0.1; } //set lower bound
		}	
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ParameterSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		interaction(size_t index_i, Real dt)
	{
		/* With Pseudo Time Step */
		this->splitting_index_[index_i] = 0;
		ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
		error_and_parameters.error_ = error_and_parameters.error_ - dt * this->residual_T_local_[index_i];
		updateStatesByError(index_i, dt, error_and_parameters);
		this->residual_k_local_[index_i] = error_and_parameters.error_;

		ErrorAndParameters<VariableType> error_and_parameters_aft = computeErrorAndParameters(index_i, dt);
		error_and_parameters_aft.error_ = error_and_parameters_aft.error_ - dt * this->residual_T_local_[index_i];
		this->residual_sp_pde_[index_i] = error_and_parameters_aft.error_;

		if (abs(this->residual_sp_pde_[index_i]) > abs(this->residual_k_local_[index_i]))
		{
			this->variable_[index_i] = this->parameter_recovery_[index_i];
			Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t& index_j = inner_neighborhood.j_[n];
				this->variable_[index_j] = this->parameter_recovery_[index_j];
			}

			ErrorAndParameters<VariableType> error_and_parameters_reverse = computeErrorAndParameters(index_i, -dt);
			error_and_parameters_reverse.error_ = error_and_parameters_reverse.error_ + dt * this->residual_T_local_[index_i];
			updateStatesByError(index_i, -dt, error_and_parameters_reverse);
			this->residual_k_local_[index_i] = error_and_parameters_reverse.error_;

			ErrorAndParameters<VariableType> error_and_parameters_aft = computeErrorAndParameters(index_i, -dt);
			error_and_parameters_aft.error_ = error_and_parameters_aft.error_ + dt * this->residual_T_local_[index_i];
			this->residual_sp_pde_[index_i] = error_and_parameters_aft.error_;

			if (abs(this->residual_sp_pde_[index_i]) > abs(this->residual_k_local_[index_i]))
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
			else if (abs(this->residual_sp_pde_[index_i]) <= abs(this->residual_k_local_[index_i]))
			{
				this->splitting_index_[index_i] = 2;
			}
		}
		else if (abs(this->residual_sp_pde_[index_i]) <= abs(this->residual_k_local_[index_i]))
		{
			this->splitting_index_[index_i] = 1;
		}
	};
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	ParameterSplittingByPDEWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                                ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		ParameterSplittingByPDEWithBoundary(ComplexRelation& complex_relation, const std::string& variable_name):
		ParameterSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, 
		                             NUM_SPECIES>(complex_relation.getInnerRelation(), variable_name),
		DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                             ContactBaseMaterialType, NUM_SPECIES>(complex_relation.getContactRelation())									 
	{
		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			boundary_heat_flux_.push_back(&(this->contact_particles_[k]->heat_flux_));
			boundary_species_.push_back(&(this->contact_particles_[k]->species_n_));
			boundary_normal_vector_.push_back(&this->contact_particles_[k]->normal_vector_);
			boundary_normal_distance_.push_back(&this->contact_particles_[k]->normal_distance_);
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> ParameterSplittingByPDEWithBoundary<BaseParticlesType, BaseMaterialType, 
		ContactBaseParticlesType, ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = ParameterSplittingByPDEInner<BaseParticlesType, 
			BaseMaterialType, VariableType, NUM_SPECIES>::computeErrorAndParameters(index_i, dt);

		VariableType& variable_i = this->variable_[index_i];
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& heat_flux_k = *(this->boundary_heat_flux_[k]);
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdLargeVec<Real>& normal_distance_k = *(this->boundary_normal_distance_[k]);
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
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	ParameterSplittingByBCInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		ParameterSplittingByBCInner(BaseInnerRelation& inner_relation, const std::string& variable_name):
		OptimizationBySplittingAlgorithmBase<BaseParticlesType, BaseMaterialType, 
		                                     VariableType, NUM_SPECIES>(inner_relation, variable_name){};*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> ParameterSplittingByBCInner<BaseParticlesType, BaseMaterialType, 
		VariableType, NUM_SPECIES>::computeErrorAndParameters(size_t index_i, Real dt)
	{
		VariableType& variable_i = this->variable_[index_i];
		ErrorAndParameters<VariableType> error_and_parameters;
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];

			VariableType variable_derivative = variable_i + this->variable_[index_j];
			Real phi_ij = this->species_modified_[index_i] - this->species_recovery_[index_j];
			Real parameter_b = -phi_ij * dW_ijV_j_ * dt;

			error_and_parameters.error_ -= variable_derivative * parameter_b;
			error_and_parameters.a_ += parameter_b;
			error_and_parameters.c_ += parameter_b * parameter_b;
		}
		return error_and_parameters;
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ParameterSplittingByBCInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters)
	{
		Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
		VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
		this->parameter_recovery_[index_i] = this->variable_[index_i];
		this->variable_[index_i] += parameter_k * error_and_parameters.a_;
		if (this->variable_[index_i] < 0.000001) { this->variable_[index_i] = 0.000001; }

		//tag the main particle modified by the BC constraint.
		this->boundary_index_[index_i] = 1;

		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];

			Real phi_ij = this->species_modified_[index_i] - this->species_recovery_[index_j];
			Real parameter_b = -phi_ij * dW_ijV_j_ * dt;

			this->parameter_recovery_[index_j] = this->variable_[index_j];
			this->variable_[index_j] += parameter_k * parameter_b;
			if (this->variable_[index_j] < 0.000001) { this->variable_[index_j] = 0.000001; }

			//tag the particles modified by the BC constraint.
			//this->boundary_index _[index_i] = 1;
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ParameterSplittingByBCInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
	     interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
			error_and_parameters.error_ = error_and_parameters.error_ + this->real_heat_flux_T_[index_i] * dt;
			this->residual_k_constrain_[index_i] = error_and_parameters.error_;
			updateStatesByError(index_i, dt, error_and_parameters);
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	ParameterSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		                               ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		ParameterSplittingByBCWithBoundary(ComplexRelation& complex_relation, const std::string& variable_name) :
		ParameterSplittingByBCInner<BaseParticlesType, BaseMaterialType, 
		                            VariableType, NUM_SPECIES>(complex_relation.getInnerRelation(), variable_name),
		DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                             ContactBaseMaterialType, NUM_SPECIES>(complex_relation.getContactRelation())
	{
		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			boundary_boundary_index_.push_back(&this->contact_particles_[k]->boundary_index_);
			boundary_heat_flux_.push_back(&(this->contact_particles_[k]->heat_flux_));
			boundary_species_.push_back(&(this->contact_particles_[k]->species_n_));
			boundary_normal_vector_.push_back(&this->contact_particles_[k]->normal_vector_);
			boundary_normal_distance_.push_back(&this->contact_particles_[k]->normal_distance_);
			boundary_parameter_recovery_.push_back(&this->contact_particles_[k]->parameter_recovery_);
			boundary_variable_.push_back(this->contact_particles_[k]->template getVariableByName<VariableType>(variable_name));
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> ParameterSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		VariableType& variable_i = this->variable_[index_i];
		ErrorAndParameters<VariableType> error_and_parameters;
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdLargeVec<Real>& normal_distance_k = *(this->boundary_normal_distance_[k]);
			StdLargeVec<Real>& heat_flux_k = *(this->boundary_heat_flux_[k]);
			StdLargeVec<VariableType>& variable_k = *(this->boundary_variable_[k]);
			StdVec<StdLargeVec<Real>>& species_k = *(boundary_species_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t& index_j = contact_neighborhood.j_[n];
				Real& dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd& e_ij_ = contact_neighborhood.e_ij_[n];

				VariableType variable_derivative = variable_k[index_j];
				Real phi_ij = this->species_modified_[index_i] - species_k[this->phi_][index_j];
				Real parameter_b = 4.0 * phi_ij * dW_ijV_j_ * this->normal_vector_[index_i].dot(e_ij_) * dt;

				error_and_parameters.error_ -= variable_derivative * parameter_b;
				error_and_parameters.c_ += parameter_b * parameter_b;
			}
		}
		return error_and_parameters;
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ParameterSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                                    ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters)
	{
		Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
		VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);

		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdLargeVec<Real>& normal_distance_k = *(this->boundary_normal_distance_[k]);
			StdLargeVec<Real>& parameter_recovery_k = *(this->boundary_parameter_recovery_[k]);
			StdLargeVec<Real>& heat_flux_k = *(this->boundary_heat_flux_[k]);
			StdLargeVec<int>& boundary_index_k = *(this->boundary_boundary_index_[k]);
			StdLargeVec<VariableType>& variable_k = *(this->boundary_variable_[k]);
			StdVec<StdLargeVec<Real>>& species_k = *(boundary_species_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t& index_j = contact_neighborhood.j_[n];
				Real& dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd& e_ij_ = contact_neighborhood.e_ij_[n];

				Real phi_ij = this->species_modified_[index_i] - species_k[this->phi_][index_j];
				Real parameter_b = 4.0 * phi_ij * dW_ijV_j_ * this->normal_vector_[index_i].dot(e_ij_) * dt;

				parameter_recovery_k[index_j] = variable_k[index_j];
				variable_k[index_j] += parameter_k * parameter_b;
				if (variable_k[index_j] < 0.000001) { variable_k[index_j] = 0.000001; }

				//tag the particles modified by the BC constraint.
				boundary_index_k[index_j] = 1;
			}
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ParameterSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		 ContactBaseMaterialType, VariableType, NUM_SPECIES>::interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
			error_and_parameters.error_ = error_and_parameters.error_ + this->real_heat_flux_T_[index_i] * dt;
			this->residual_k_constrain_[index_i] = error_and_parameters.error_;
			updateStatesByError(index_i, dt, error_and_parameters);
		}
	}*/
	//=================================================================================================//
	/*template <typename ParameterSplittingType, typename BaseBodyRelationType, typename VariableType>
	UpdateParameterBCResidual<ParameterSplittingType, BaseBodyRelationType, VariableType>::
		UpdateParameterBCResidual(BaseBodyRelationType& body_relation, const std::string& variable_name) :
		ParameterSplittingType(body_relation, variable_name) {}*/
	//=================================================================================================//
	/*template <typename ParameterSplittingType, typename BaseBodyRelationType, typename VariableType>
	void UpdateParameterBCResidual<ParameterSplittingType, BaseBodyRelationType, VariableType>::
		interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			ErrorAndParameters<VariableType> error_and_parameters = this->computeErrorAndParameters(index_i, dt);
			this->real_heat_flux_k_[index_i] = -error_and_parameters.error_ / dt;
		}
	}*/
	//=================================================================================================//
}
#endif DIFFUSION_SPLITTING_PARAMETER_HPP
