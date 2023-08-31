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
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	TemperatureSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		TemperatureSplittingByPDEInner(BaseInnerRelation& inner_relation, const std::string& variable_name) :
		OptimizationBySplittingAlgorithmBase<BaseParticlesType, BaseMaterialType, 
		                                     VariableType, NUM_SPECIES>(inner_relation, variable_name) {};
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> TemperatureSplittingByPDEInner<BaseParticlesType, BaseMaterialType, 
		VariableType, NUM_SPECIES>::computeErrorAndParameters(size_t index_i, Real dt)
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
			Real diff_coff_ij = this->species_diffusion_[this->phi_]->getInterParticleDiffusionCoff(index_i, index_j, e_ij_);
			Real parameter_b = 2.0 * diff_coff_ij * inner_neighborhood.dW_ijV_j_[n] * dt / r_ij_;

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
	void TemperatureSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
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
	
			Real diff_coff_ij = this->species_diffusion_[this->phi_]->getInterParticleDiffusionCoff(index_i, index_j, e_ij_);
			Real parameter_b = 2.0 * diff_coff_ij * inner_neighborhood.dW_ijV_j_[n] * dt / r_ij_;
			this->variable_[index_j] -= parameter_k * parameter_b;
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void TemperatureSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		interaction(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
		updateStatesByError(index_i, dt, error_and_parameters);
		this->residual_T_local_[index_i] = error_and_parameters.error_;
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
			  class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	TemperatureSplittingByPDEWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                                  ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		TemperatureSplittingByPDEWithBoundary(ComplexRelation & complex_relation, const std::string & variable_name):
		TemperatureSplittingByPDEInner<BaseParticlesType, BaseMaterialType, VariableType, 
		                               NUM_SPECIES>(complex_relation.getInnerRelation(), variable_name),
		DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                             ContactBaseMaterialType, NUM_SPECIES>(complex_relation.getContactRelation())
	{
		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			boundary_normal_vector_.push_back(&this->contact_particles_[k]->normal_vector_);
			boundary_heat_flux_.push_back(&this->contact_particles_[k]->heat_flux_);
			boundary_normal_distance_.push_back(&(this->contact_particles_[k]->normal_distance_));
			boundary_variable_.push_back(this->contact_particles_[k]->template getVariableByName<VariableType>(variable_name));
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> TemperatureSplittingByPDEWithBoundary<BaseParticlesType, BaseMaterialType, 
		ContactBaseParticlesType, ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		ErrorAndParameters<VariableType> error_and_parameters = TemperatureSplittingByPDEInner<BaseParticlesType,
			BaseMaterialType, VariableType, NUM_SPECIES>::computeErrorAndParameters(index_i, dt);

		VariableType& variable_i = this->variable_[index_i];
		/* contact interaction. */
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& heat_flux_k = *(this->boundary_heat_flux_[k]);
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdLargeVec<Real>& normal_distance_k = *(this->boundary_normal_distance_[k]);
			StdLargeVec<VariableType>& variable_k = *(this->boundary_variable_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t& index_j = contact_neighborhood.j_[n];
				
				if (variable_k[index_j] > 0.0)
				{
					// linear projection
					VariableType variable_derivative = (variable_i - variable_k[index_j]);
					Real diff_coff_ij = this->species_diffusion_[this->phi_]->getDiffusionCoffWithBoundary(index_i);
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
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	TemperatureSplittingByBCInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		TemperatureSplittingByBCInner(BaseInnerRelation& inner_relation, const std::string& variable_name) :
		OptimizationBySplittingAlgorithmBase<BaseParticlesType, BaseMaterialType,
		                                     VariableType, NUM_SPECIES>(inner_relation, variable_name) {};*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> TemperatureSplittingByBCInner<BaseParticlesType, BaseMaterialType,
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
			Vecd& e_ij_ = inner_neighborhood.e_ij_[n];

			VariableType variable_derivative = (variable_i - this->variable_[index_j]);
			Real diff_coff_ij = this->species_diffusion_[this->phi_]->getInterParticleDiffusionCoff(index_i, index_j, e_ij_);
			Real parameter_b = -2.0 * diff_coff_ij * dW_ijV_j_ * dt;

			error_and_parameters.error_ -= variable_derivative * parameter_b;
			error_and_parameters.a_ += parameter_b;
			error_and_parameters.c_ += parameter_b * parameter_b;
		}
		error_and_parameters.error_ = error_and_parameters.error_ + this->heat_flux_[index_i] * dt;
		return error_and_parameters;
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void TemperatureSplittingByBCInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters)
	{
		Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
		VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
		this->variable_[index_i] += parameter_k * error_and_parameters.a_;

		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t& index_j = inner_neighborhood.j_[n];
			Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real& r_ij_ = inner_neighborhood.r_ij_[n];
			Vecd& e_ij_ = inner_neighborhood.e_ij_[n];

			Real diff_coff_ij = this->species_diffusion_[this->phi_]->getInterParticleDiffusionCoff(index_i, index_j, e_ij_);
			Real parameter_b = -2.0 * diff_coff_ij * dW_ijV_j_ * dt;
			this->variable_[index_j] -= parameter_k * parameter_b;
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void TemperatureSplittingByBCInner<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
			updateStatesByError(index_i, dt, error_and_parameters);
			this->residual_T_constrain_[index_i] = error_and_parameters.error_;
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
		TemperatureSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		                                     ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		TemperatureSplittingByBCWithBoundary(ComplexRelation& complex_relation, const std::string& variable_name) :
		TemperatureSplittingByBCInner<BaseParticlesType, BaseMaterialType,
		                              VariableType, NUM_SPECIES>(complex_relation.getInnerRelation(), variable_name),
		DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		ContactBaseMaterialType, NUM_SPECIES>(complex_relation.getContactRelation())
	{
		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			boundary_normal_vector_.push_back(&this->contact_particles_[k]->normal_vector_);
			boundary_normal_distance_.push_back(&this->contact_particles_[k]->normal_distance_);
			boundary_variable_.push_back(this->contact_particles_[k]->template getVariableByName<VariableType>(variable_name));
			boundary_diffusivity_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("ThermalDiffusivity"));
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	ErrorAndParameters<VariableType> TemperatureSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		computeErrorAndParameters(size_t index_i, Real dt)
	{
		VariableType& variable_i = this->variable_[index_i];
		ErrorAndParameters<VariableType> error_and_parameters;

		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& diffusivity_k = *(this->boundary_diffusivity_[k]);
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdLargeVec<Real>& normal_distance_k = *(this->boundary_normal_distance_[k]);
			StdLargeVec<VariableType>& variable_k = *(this->boundary_variable_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t& index_j = contact_neighborhood.j_[n];
				Real& dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Real& r_ij_ = contact_neighborhood.r_ij_[n];
				Vecd& e_ij_ = contact_neighborhood.e_ij_[n];

				VariableType variable_derivative = (variable_i - variable_k[index_j]);
				Real diff_coff_ij = diffusivity_k[index_j];
				Real parameter_b = 4.0 * diff_coff_ij * dW_ijV_j_ * this->normal_vector_[index_i].dot(e_ij_) * dt;

				error_and_parameters.error_ -= variable_derivative * parameter_b;
				error_and_parameters.a_ += parameter_b;
				error_and_parameters.c_ += parameter_b * parameter_b;
			}
		}
		error_and_parameters.error_ = error_and_parameters.error_ + this->heat_flux_[index_i] * dt;
		return error_and_parameters;
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	void TemperatureSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		                                      ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		updateStatesByError(size_t index_i, Real dt, const ErrorAndParameters<VariableType>& error_and_parameters)
	{
		Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
		VariableType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
		this->variable_[index_i] += parameter_k * error_and_parameters.a_;

		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real>& diffusivity_k = *(this->boundary_diffusivity_[k]);
			StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
			StdLargeVec<Real>& normal_distance_k = *(this->boundary_normal_distance_[k]);
			StdLargeVec<VariableType>& variable_k = *(this->boundary_variable_[k]);

			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t& index_j = contact_neighborhood.j_[n];
				Real& dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Real& r_ij_ = contact_neighborhood.r_ij_[n];
				Vecd& e_ij_ = contact_neighborhood.e_ij_[n];

				Real diff_coff_ij = diffusivity_k[index_j];
				Real parameter_b = 4.0 * diff_coff_ij * dW_ijV_j_ * this->normal_vector_[index_i].dot(e_ij_) * dt;
				variable_k[index_j] -= parameter_k * parameter_b;
			}
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		      class ContactBaseMaterialType, typename VariableType, int NUM_SPECIES>
	void TemperatureSplittingByBCWithBoundary<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		ContactBaseMaterialType, VariableType, NUM_SPECIES>::interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			ErrorAndParameters<VariableType> error_and_parameters = computeErrorAndParameters(index_i, dt);
			updateStatesByError(index_i, dt, error_and_parameters);
			this->residual_T_constrain_[index_i] = error_and_parameters.error_;
		}
	}*/
	//=================================================================================================//
	/*template <typename TemperatureConstrainType, typename BaseBodyRelationType, typename VariableType>
	UpdateTemperatureBCResidual<TemperatureConstrainType, BaseBodyRelationType, VariableType>::
		UpdateTemperatureBCResidual(BaseBodyRelationType& body_relation, const std::string& variable_name) :
		TemperatureConstrainType(body_relation, variable_name) {};*/
	//=================================================================================================//
	/*template <typename TemperatureConstrainType, typename BaseBodyRelationType, typename VariableType>
	void  UpdateTemperatureBCResidual<TemperatureConstrainType, BaseBodyRelationType, VariableType>::
		interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			ErrorAndParameters<VariableType> error_and_parameters = this->computeErrorAndParameters(index_i, dt);
			this->real_heat_flux_T_[index_i] = this->heat_flux_[index_i] - error_and_parameters.error_ / dt;
		}
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	UpdateBoundaryParticleTemperature<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		UpdateBoundaryParticleTemperature(BaseInnerRelation& inner_relation, const std::string& variable_name) :
		LocalDynamics(inner_relation.sph_body_),
		DiffusionReactionInnerData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(inner_relation),
		sigma0_(sph_body_.sph_adaptation_->ReferenceNumberDensity()),
		heat_flux_(this->particles_->heat_flux_), normal_vector_(this->particles_->normal_vector_),
		variable_(*this->particles_->template getVariableByName<VariableType>(variable_name))
	{
		phi_ = this->particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}*/
	//=================================================================================================//
	/*template <class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void UpdateBoundaryParticleTemperature<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			Real parameter_i, parameter_j;
			VariableType& variable_i = this->variable_[index_i];
			Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];

			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t& index_j = inner_neighborhood.j_[n];
				Real& dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
				Real& r_ij_ = inner_neighborhood.r_ij_[n];
				Vecd& e_ij_ = inner_neighborhood.e_ij_[n];

				//Here the temperature gradient involves all neighboring particles, not judged from the normal vector,
				  and the coefficient here needs to be adjusted according to the number of its neighboring.
				Real diff_coff_ij = this->species_diffusion_[this->phi_]->getDiffusionCoffWithBoundary(index_i);
				Real parameter_b = -2.0 * diff_coff_ij * dW_ijV_j_;

				parameter_i += parameter_b;
				parameter_j += variable_[index_j] * parameter_b;
			}
			variable_[index_i] = (heat_flux_[index_i] + parameter_j) / (parameter_i + TinyReal);
		}
	}*/
	//=================================================================================================//
	/*template<class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		     class ContactBaseMaterialType, class VariableType, int NUM_SPECIES>
	UpdateBoundaryParticleTemperatureDummy<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                                   ContactBaseMaterialType, VariableType, NUM_SPECIES>::
		UpdateBoundaryParticleTemperatureDummy(ComplexRelation& body_complex_relation, const std::string& variable_name):
		UpdateBoundaryParticleTemperature<BaseParticlesType, BaseMaterialType, VariableType,  
		                                  NUM_SPECIES>(body_complex_relation.getInnerRelation(), variable_name),
		DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                             ContactBaseMaterialType, NUM_SPECIES>(body_complex_relation.getContactRelation())
	{
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			boundary_normal_vector_.push_back(&this->contact_particles_[k]->normal_vector_);
			boundary_variable_.push_back(this->contact_particles_[k]->template getVariableByName<VariableType>(variable_name));
			boundary_diffusivity_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("ThermalDiffusivity"));
		}
	}*/
	//=================================================================================================//
	/*template<class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType,
		     class ContactBaseMaterialType, class VariableType, int NUM_SPECIES>
    void UpdateBoundaryParticleTemperatureDummy<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		ContactBaseMaterialType, VariableType, NUM_SPECIES>::interaction(size_t index_i, Real dt)
	{
		if (this->heat_flux_[index_i] != 0)
		{
			Real parameter_i, parameter_j, W;
			for (size_t k = 0; k != this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& normal_vector_k = *(this->boundary_normal_vector_[k]);
				StdLargeVec<Real>& variable_k = *(this->boundary_variable_[k]);
				StdLargeVec<Real>& diffusivity_k = *(this->boundary_diffusivity_[k]);
				Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];

				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t& index_j = contact_neighborhood.j_[n];
					Real& dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
					Real& r_ij_ = contact_neighborhood.r_ij_[n];
					Vecd& e_ij_ = contact_neighborhood.e_ij_[n];

					//Here the temperature gradient only involves the thermal domain particles, judged by the 
					  normal vector, and the coefficient here needs to be adjusted based on the correction.
					Real parameter_b = 4.0 * diffusivity_k[index_j] * dW_ijV_j_ * this->normal_vector_[index_i].dot(e_ij_);
					parameter_i += parameter_b;
					parameter_j += variable_k[index_j] * parameter_b;
				}
			}
			this->variable_[index_i] = (this->heat_flux_[index_i] + parameter_j) / (parameter_i + TinyReal);
		}
	}*/
	//=================================================================================================//
}

#endif DIFFUSION_SPLITTING_STATE_HPP