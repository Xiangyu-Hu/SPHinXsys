/**
* @file 	diffusion_optimization_common.hpp
* @author	Bo Zhang and Xiangyu Hu
*/

#ifndef DIFFUSION_OPTIMIZATION_COMMON_HPP
#define DIFFUSION_OPTIMIZATION_COMMON_HPP

#include "diffusion_optimization_common.h"

namespace SPH
{
	//=================================================================================================//
	template<class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType, 
		     class ContactBaseMaterialType, int NUM_SPECIES>
	UpdateUnitNormalVector<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType,
		                   ContactBaseMaterialType, NUM_SPECIES>::
		UpdateUnitNormalVector(ComplexRelation& body_complex_relation) :
		LocalDynamics(body_complex_relation.getInnerRelation().sph_body_),
		DiffusionReactionInnerData<BaseParticlesType, BaseMaterialType, 
		                           NUM_SPECIES>(body_complex_relation.getInnerRelation()),
		DiffusionReactionContactData<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, ContactBaseMaterialType, 
		                             NUM_SPECIES>(body_complex_relation.getContactRelation()),
		normal_vector_(this->particles_->normal_vector_) {}
	//=================================================================================================//
	template<class BaseParticlesType, class BaseMaterialType, class ContactBaseParticlesType, 
		     class ContactBaseMaterialType, int NUM_SPECIES>
	void UpdateUnitNormalVector<BaseParticlesType, BaseMaterialType, ContactBaseParticlesType, 
		                        ContactBaseMaterialType, NUM_SPECIES>::interaction(size_t index_i, Real dt)
	{
		for (size_t k = 0; k != this->contact_configuration_.size(); ++k)
		{
			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				Real &dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd& e_ij_ = contact_neighborhood.e_ij_[n];
				normal_vector_[index_i] += dW_ijV_j_ * e_ij_;
			}
		}
		normal_vector_[index_i] = normal_vector_[index_i] / (normal_vector_[index_i].norm() + TinyReal);
	};
	//=================================================================================================//
	template<class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	UpdateNormalDistance<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		UpdateNormalDistance(BaseInnerRelation& inner_relation) :
		LocalDynamics(inner_relation.sph_body_),
		DiffusionReactionInnerData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(inner_relation),
		W0_(sph_body_.sph_adaptation_->getKernel()->W0(zero_vec)), 
		sigma0_(sph_body_.sph_adaptation_->ReferenceNumberDensity()),
		cutoff_radius_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
		normal_distance_(this->particles_->normal_distance_) {};
	//=================================================================================================//
	template<class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void UpdateNormalDistance<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		 interaction(size_t index_i, Real dt)
	{
		Real sigma_(W0_);
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			Real& W_ij_ = inner_neighborhood.W_ij_[n];
			sigma_ += W_ij_;
		}
		normal_distance_[index_i] = cutoff_radius_ * (2 * sigma_ / sigma0_ - 1);
	}
	//=================================================================================================//
	template<class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	ThermalDiffusivityConstrain<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		ThermalDiffusivityConstrain(SPHBody& diffusion_body, const std::string& variable_name,
			                        Real initial_thermal_diffusivity):
		LocalDynamics(diffusion_body), 
		DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(diffusion_body),
		initial_thermal_diffusivity_(initial_thermal_diffusivity), new_averaged_thermal_diffusivity_(0.0),
		local_thermal_conductivity_(*this->particles_->template getVariableByName<VariableType>(variable_name)) {};
	//=================================================================================================//
	template<class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ThermalDiffusivityConstrain<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		UpdateAveragedParameter(Real new_averaged_thermal_diffusivity)
	{
		new_averaged_thermal_diffusivity_ = new_averaged_thermal_diffusivity;
	}
	//=================================================================================================//
	template<class BaseParticlesType, class BaseMaterialType, typename VariableType, int NUM_SPECIES>
	void ThermalDiffusivityConstrain<BaseParticlesType, BaseMaterialType, VariableType, NUM_SPECIES>::
		update(size_t index_i, Real dt)
	{
		local_thermal_conductivity_[index_i] = local_thermal_conductivity_[index_i] *
			initial_thermal_diffusivity_ / new_averaged_thermal_diffusivity_;
	}
	//=================================================================================================//
}
#endif DIFFUSION_OPTIMIZATION_COMMON_HPP