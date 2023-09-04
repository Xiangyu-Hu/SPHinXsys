/**
* @file 	diffusion_splitting_state.hpp
* @author	Bo Zhang and Xiangyu Hu
*/

#ifndef DIFFUSION_SPLITTING_STATE_HPP
#define DIFFUSION_SPLITTING_STATE_HPP

#include "diffusion_splitting_state.h"

namespace SPH
{
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	TemperatureSplittingByPDEInner<ParticlesType, VariableType>::
		TemperatureSplittingByPDEInner(BaseInnerRelation &inner_relation, const std::string &variable_name) : 
		OptimizationBySplittingAlgorithmBase<ParticlesType, VariableType>(inner_relation, variable_name){};
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
    ErrorAndParameters<VariableType> TemperatureSplittingByPDEInner<ParticlesType, VariableType>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		VariableType& variable_i = this->variable_[index_i];
		ErrorAndParameters<VariableType> error_and_parameters;
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& r_ij_ = inner_neighborhood.r_ij_[n];
			Vecd& e_ij_ = inner_neighborhood.e_ij_[n];

			// linear projection
			VariableType variable_derivative = (variable_i - this->variable_[index_j]);
			Real diff_coff_ij = this->all_diffusion_[this->phi_]->getInterParticleDiffusionCoeff(index_i, index_j, e_ij_);
			Real parameter_b = 2.0 * diff_coff_ij * inner_neighborhood.dW_ijV_j_[n] * dt / r_ij_;

			error_and_parameters.error_ -= variable_derivative * parameter_b;
			error_and_parameters.a_ += parameter_b;
			error_and_parameters.c_ += parameter_b * parameter_b;
		}
		error_and_parameters.a_ -= 1;
		error_and_parameters.error_ -= this->heat_source_[index_i] * dt;
		return error_and_parameters;
	};
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	void TemperatureSplittingByPDEInner<ParticlesType, VariableType>::
		updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters)
	{
		Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
		VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
		this->variable_[index_i] += parameter_k * error_and_parameters.a_;

		Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& r_ij_ = inner_neighborhood.r_ij_[n];
			Vecd& e_ij_ = inner_neighborhood.e_ij_[n];
	
			Real diff_coff_ij = this->all_diffusion_[this->phi_]->getInterParticleDiffusionCoeff(index_i, index_j, e_ij_);
			Real parameter_b = 2.0 * diff_coff_ij * inner_neighborhood.dW_ijV_j_[n] * dt / r_ij_;
			this->variable_[index_j] -= parameter_k * parameter_b;
		}
	}
	//=================================================================================================//
	template <class ParticlesType, typename VariableType>
	void TemperatureSplittingByPDEInner<ParticlesType, VariableType>::
		interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
		updateStatesByError(index_i, dt, error_and_parameters);
		this->residual_T_local_[index_i] = error_and_parameters.error_;
	}
	//=================================================================================================//
	template <class ParticlesType, class ContactParticlesType, typename VariableType>
	TemperatureSplittingByPDEWithBoundary<ParticlesType, ContactParticlesType, VariableType>::
		TemperatureSplittingByPDEWithBoundary(ComplexRelation &complex_relation, const std::string &variable_name):
		TemperatureSplittingByPDEInner<ParticlesType, VariableType>(complex_relation.getInnerRelation(), variable_name),
		DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>(complex_relation.getContactRelation())
    {
		boundary_heat_flux_.resize(this->contact_particles_.size());
		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			boundary_normal_vector_.push_back(&this->contact_particles_[k]->n_);			 
			boundary_variable_.push_back(this->contact_particles_[k]->template getVariableByName<VariableType>(variable_name));
			boundary_heat_flux_[k] = this->contact_particles_[k]->template registerSharedVariable<Real>("HeatFlux");
		}
	};
	//=================================================================================================//
	template <class ParticlesType, class ContactParticlesType, typename VariableType>
	ErrorAndParameters<VariableType> TemperatureSplittingByPDEWithBoundary<ParticlesType, ContactParticlesType, VariableType>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = TemperatureSplittingByPDEInner<ParticlesType, VariableType>
			::computeErrorAndParameters(index_i, dt);

		VariableType& variable_i = this->variable_[index_i];
		/* contact interaction. */
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& heat_flux_k = *(this->boundary_heat_flux_[k]);
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdLargeVec<VariableType>& variable_k = *(this->boundary_variable_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t& index_j = contact_neighborhood.j_[n];
				
				if (variable_k[index_j] > 0.0)
				{
					// linear projection
					VariableType variable_derivative = (variable_i - variable_k[index_j]);
					Real diff_coff_ij = this->all_diffusion_[this->phi_]->getDiffusionCoeffWithBoundary(index_i);
					Real parameter_b = 2.0 * diff_coff_ij * contact_neighborhood.dW_ijV_j_[n] * dt / contact_neighborhood.r_ij_[n];

					error_and_parameters.error_ -= variable_derivative * parameter_b;
					error_and_parameters.a_ += parameter_b;
				}

				Vecd n_ij = this->normal_vector_[index_i] - normal_vector_k[index_j];
				error_and_parameters.error_ -= heat_flux_k[index_j] * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n].dot(n_ij) * dt;
			}
		}
		return error_and_parameters;
	}
	//=================================================================================================//
	template <typename TemperatureSplittingType, typename BaseBodyRelationType, typename VariableType>
	UpdateTemperaturePDEResidual<TemperatureSplittingType, BaseBodyRelationType, VariableType>::
		UpdateTemperaturePDEResidual(BaseBodyRelationType& body_relation, const std::string& variable_name) :
		TemperatureSplittingType(body_relation, variable_name) {};
	//=================================================================================================//
	template <typename TemperatureSplittingType, typename BaseBodyRelationType, typename VariableType>
	void UpdateTemperaturePDEResidual<TemperatureSplittingType, BaseBodyRelationType, VariableType>::
		interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = this->computeErrorAndParameters(index_i, dt);
		this->residual_T_global_[index_i] = error_and_parameters.error_;
	}
	//=================================================================================================//
}

#endif //DIFFUSION_SPLITTING_STATE_HPP