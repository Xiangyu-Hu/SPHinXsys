/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	particle_dynamics_diffusion_reaction.hpp
 * @brief 	This is the particle dynamics applicable for all type bodies
 * 			TODO: there is an issue on applying corrected configuration for contact bodies.. 
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef  EXTRA_TWO_PHASE_HEAT_TRANSFER_HPP
#define  EXTRA_TWO_PHASE_HEAT_TRANSFER_HPP

#include "extra_two_phase_heat_transfer.h"

namespace SPH
{
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	TwoPhaseDiffusionReactionInitialCondition<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		TwoPhaseDiffusionReactionInitialCondition(SPHBody &sph_body)
		: LocalDynamics(sph_body),
		TwoPhaseDiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(sph_body),
		  pos_(this->particles_->pos_), species_n_(this->particles_->species_n_),
		  thermal_conductivity_(this->particles_->thermal_conductivity_) {}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	TwoPhaseGetDiffusionTimeStepSize<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		TwoPhaseGetDiffusionTimeStepSize(SPHBody &sph_body)
		: BaseDynamics<Real>(),
		TwoPhaseDiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(sph_body)
	{
		Real smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
		diff_time_step_ = this->particles_->diffusion_reaction_material_.getDiffusionTimeStepSize(smoothing_length);
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		TwoPhaseRelaxationOfAllDiffusionSpeciesInner(BaseInnerRelation &inner_relation)
		: LocalDynamics(inner_relation.sph_body_),
		TwoPhaseDiffusionReactionInnerData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(inner_relation),
		  diffusion_reaction_material_(this->particles_->diffusion_reaction_material_),
		  species_n_(this->particles_->species_n_),
		  diffusion_dt_(this->particles_->diffusion_dt_)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		initializeDiffusionChangeRate(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			diffusion_dt_[m][particle_i] = 0;
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			Real diff_coff_ij = species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
			size_t l = species_diffusion_[m]->gradient_species_index_;
			Real phi_ij = species_n_[l][particle_i] - species_n_[l][particle_j];
			diffusion_dt_[m][particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_n_[k][particle_i] += dt * diffusion_dt_[m][particle_i];
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		interaction(size_t index_i, Real dt)
	{
		TwoPhaseDiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES> *particles = this->particles_;
		Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];

		initializeDiffusionChangeRate(index_i);
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real r_ij_ = inner_neighborhood.r_ij_[n];
			Vecd &e_ij = inner_neighborhood.e_ij_[n];

			const Vecd &grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
			Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
			getDiffusionChangeRate(index_i, index_j, e_ij, area_ij);
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		update(size_t index_i, Real dt)
	{
		updateSpeciesDiffusion(index_i, dt);
	}
	
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	TwoPhaseInitializationRK<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		TwoPhaseInitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &species_s)
		: LocalDynamics(sph_body),
		TwoPhaseDiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(sph_body),
		  species_n_(this->particles_->species_n_), species_s_(species_s)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void TwoPhaseInitializationRK<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		initializeIntermediateValue(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_s_[m][particle_i] = species_n_[k][particle_i];
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType, int NUM_SPECIES>
	void TwoPhaseInitializationRK<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::
		update(size_t index_i, Real dt)
	{
		initializeIntermediateValue(index_i);
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
		class ContactBaseParticlesType, class ContactBaseMaterialType, int NUM_SPECIES>
	TwoPhaseRelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>::
		TwoPhaseRelaxationOfAllDiffusionSpeciesComplex(ComplexRelation& complex_relation)
		: TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>(complex_relation.getInnerRelation()),
		TwoPhaseDiffusionReactionContactData<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>(complex_relation.getContactRelation()),
		species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_),
		external_diffusion_dt_(this->particles_->external_diffusion_dt_), thermal_conductivity_(this->particles_->thermal_conductivity_),
		external_diffusion_dt_sum_(this->particles_->external_diffusion_dt_sum_)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();

		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			contact_species_n_.push_back(&(this->contact_particles_[k]->species_n_));
			contact_thermal_conductivity_.push_back(&(this->contact_particles_[k]->thermal_conductivity_));
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
		class ContactBaseParticlesType, class ContactBaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>::
		initializeExternalDiffusionChangeRate(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			external_diffusion_dt_[particle_i] = 0;
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
		class ContactBaseParticlesType, class ContactBaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>::
		getDiffusionChangeRateContact(size_t particle_i, size_t particle_j, Vecd& e_ij,
			Real surface_area_ij, const StdVec<StdLargeVec<Real>>& species_n_k, const StdLargeVec<Real>& thermal_conductivity_k)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			Real diff_coff_ij = 2.0 * species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij)
				* thermal_conductivity_k[particle_j] / (thermal_conductivity_[particle_i] + thermal_conductivity_k[particle_j]);
			//Real diff_coff_ij = species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
			size_t l = species_diffusion_[m]->gradient_species_index_;
			Real phi_ij = species_n_[l][particle_i] - species_n_k[l][particle_j];
			external_diffusion_dt_[particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
			diffusion_dt_[m][particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
		class ContactBaseParticlesType, class ContactBaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_n_[k][particle_i] += dt * diffusion_dt_[m][particle_i];
			external_diffusion_dt_sum_[particle_i] = external_diffusion_dt_[particle_i];
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
		class ContactBaseParticlesType, class ContactBaseMaterialType, int NUM_SPECIES>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
		ContactBaseParticlesType, ContactBaseMaterialType, NUM_SPECIES>::
		interaction(size_t index_i, Real dt)
	{
		TwoPhaseRelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType, NUM_SPECIES>::interaction(index_i, dt);
		TwoPhaseDiffusionReactionParticles<BaseParticlesType, BaseMaterialType, NUM_SPECIES>* particles = this->particles_;
		initializeExternalDiffusionChangeRate(index_i);
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdVec<StdLargeVec<Real>>& species_n_k = *(contact_species_n_[k]);
			StdLargeVec<Real>& thermal_conductivity_k = *(contact_thermal_conductivity_[k]);
			Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real r_ij_ = contact_neighborhood.r_ij_[n];
				Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd& e_ij = contact_neighborhood.e_ij_[n];

				const Vecd& grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
				Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
				getDiffusionChangeRateContact(index_i, index_j, e_ij, area_ij, species_n_k, thermal_conductivity_k);

			}
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	TwoPhaseSecondStageRK2<FirstStageType>::
		TwoPhaseSecondStageRK2(typename FirstStageType::BodyRelationType& body_relation,
			StdVec<StdLargeVec<Real>>& species_s)
		: FirstStageType(body_relation),
		species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_),
		species_s_(species_s)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class FirstStageType>
	void TwoPhaseSecondStageRK2<FirstStageType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < this->species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_n_[k][particle_i] = 0.5 * species_s_[m][particle_i] +
				0.5 * (species_n_[k][particle_i] + dt * diffusion_dt_[m][particle_i]);
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	TwoPhaseRelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::
		TwoPhaseRelaxationOfAllDiffusionSpeciesRK2(typename FirstStageType::BodyRelationType& body_relation)
		: BaseDynamics<void>(), rk2_initialization_(body_relation.sph_body_, species_s_),
		rk2_1st_stage_(body_relation), rk2_2nd_stage_(body_relation, species_s_)
	{
		StdVec<BaseDiffusion*> species_diffusion_ = rk2_1st_stage_.diffusion_reaction_material_.SpeciesDiffusion();

		size_t number_of_diffusion_species = species_diffusion_.size();
		species_s_.resize(number_of_diffusion_species);
		for (size_t m = 0; m < number_of_diffusion_species; ++m)
		{
			// the size should be the same as that in the base particles
			species_s_[m].resize(rk2_1st_stage_.getParticles()->real_particles_bound_);
			// register data in base particles
			std::get<0>(rk2_1st_stage_.getParticles()->all_particle_data_).push_back(&species_s_[m]);
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::exec(Real dt)
	{
		rk2_initialization_.exec();
		rk2_1st_stage_.exec(dt);
		rk2_2nd_stage_.exec(dt);
	}
	//=================================================================================================//
	template <class FirstStageType>
	void TwoPhaseRelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::parallel_exec(Real dt)
	{
		rk2_initialization_.parallel_exec();
		rk2_1st_stage_.parallel_exec(dt);
		rk2_2nd_stage_.parallel_exec(dt);
	}
	//=================================================================================================//
}
#endif // EXTRA_TWO_PHASE_HEAT_TRANSFER_HPP