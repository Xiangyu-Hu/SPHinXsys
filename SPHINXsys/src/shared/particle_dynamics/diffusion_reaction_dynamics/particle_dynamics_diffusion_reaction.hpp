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

#ifndef PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP
#define PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP

#include "particle_dynamics_diffusion_reaction.h"

namespace SPH
{
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	DiffusionReactionInitialCondition<DiffusionReactionParticlesType>::
		DiffusionReactionInitialCondition(SPHBody &sph_body)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<DiffusionReactionParticlesType>(sph_body),
		  pos_(this->particles_->pos_), all_species_(this->particles_->all_species_) {}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	GetDiffusionTimeStepSize<DiffusionReactionParticlesType>::
		GetDiffusionTimeStepSize(SPHBody &sph_body)
		: BaseDynamics<Real>(sph_body),
		  DiffusionReactionSimpleData<DiffusionReactionParticlesType>(sph_body)
	{
		Real smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
		diff_time_step_ = this->particles_->diffusion_reaction_material_
							  .getDiffusionTimeStepSize(smoothing_length);
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		RelaxationOfAllDiffusionSpeciesInner(BaseInnerRelation &inner_relation)
		: LocalDynamics(inner_relation.getSPHBody()),
		  DiffusionReactionInnerData<DiffusionReactionParticlesType>(inner_relation),
		  diffusion_reaction_material_(this->particles_->diffusion_reaction_material_),
		  all_species_(this->particles_->all_species_),
		  diffusion_dt_(this->particles_->diffusion_dt_)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		initializeDiffusionChangeRate(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			diffusion_dt_[m][particle_i] = 0;
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			Real diff_coff_ij = species_diffusion_[m]
									->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
			size_t l = species_diffusion_[m]->gradient_species_index_;
			Real phi_ij = all_species_[l][particle_i] - all_species_[l][particle_j];
			diffusion_dt_[m][particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			all_species_[k][particle_i] += dt * diffusion_dt_[m][particle_i];
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		interaction(size_t index_i, Real dt)
	{
		DiffusionReactionParticlesType *particles = this->particles_;
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
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		update(size_t index_i, Real dt)
	{
		updateSpeciesDiffusion(index_i, dt);
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	RelaxationOfAllDiffusionSpeciesComplex<
		DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		RelaxationOfAllDiffusionSpeciesComplex(ComplexRelation &complex_relation)
		: RelaxationOfAllDiffusionSpeciesInner<
			  DiffusionReactionParticlesType>(complex_relation.getInnerRelation()),
		  DiffusionReactionContactData<
			  DiffusionReactionParticlesType,
			  ContactDiffusionReactionParticlesType>(complex_relation.getContactRelation()),
		  all_species_(this->particles_->all_species_), diffusion_dt_(this->particles_->diffusion_dt_)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();

		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			contact_species_n_.push_back(&(this->contact_particles_[k]->all_species_));
		}

		StdVec<std::string> &all_species_names = this->particles_->AllSpeciesNames();
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t l = species_diffusion_[m]->gradient_species_index_;
			const std::string &inner_species = all_species_names[l];
			for (size_t k = 0; k != this->contact_particles_.size(); ++k)
			{
				StdVec<std::string> &all_species_names_k =
					this->contact_particles_[k]->AllSpeciesNames();
				const std::string &contact_species_k = all_species_names_k[l];
				if (inner_species != contact_species_k)
				{
					std::cout << "\n Error: inner species '" << inner_species
							  << "' and contact species '" << contact_species_k << "' not match! " << std::endl;
					std::cout << __FILE__ << ':' << __LINE__ << std::endl;
					exit(1);
				}
			}
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesComplex<
		DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		getDiffusionChangeRateContact(size_t particle_i, size_t particle_j, Vecd &e_ij,
									  Real surface_area_ij, const StdVec<StdLargeVec<Real>> &species_n_k)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			Real diff_coff_ij = species_diffusion_[m]
									->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
			size_t l = species_diffusion_[m]->gradient_species_index_;
			Real phi_ij = all_species_[l][particle_i] - species_n_k[l][particle_j];
			diffusion_dt_[m][particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesComplex<
		DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		interaction(size_t index_i, Real dt)
	{
		RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::interaction(index_i, dt);
		DiffusionReactionParticlesType *particles = this->particles_;

		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdVec<StdLargeVec<Real>> &species_n_k = *(contact_species_n_[k]);

			Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real r_ij_ = contact_neighborhood.r_ij_[n];
				Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd &e_ij = contact_neighborhood.e_ij_[n];

				const Vecd &grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
				Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
				getDiffusionChangeRateContact(index_i, index_j, e_ij, area_ij, species_n_k);
			}
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	InitializationRK<DiffusionReactionParticlesType>::
		InitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &species_s)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<DiffusionReactionParticlesType>(sph_body),
		  all_species_(this->particles_->all_species_), species_s_(species_s)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void InitializationRK<DiffusionReactionParticlesType>::
		initializeIntermediateValue(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_s_[m][particle_i] = all_species_[k][particle_i];
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void InitializationRK<DiffusionReactionParticlesType>::
		update(size_t index_i, Real dt)
	{
		initializeIntermediateValue(index_i);
	}
	//=================================================================================================//
	template <class FirstStageType>
	SecondStageRK2<FirstStageType>::
		SecondStageRK2(typename FirstStageType::BodyRelationType &body_relation,
					   StdVec<StdLargeVec<Real>> &species_s)
		: FirstStageType(body_relation),
		  all_species_(this->particles_->all_species_), diffusion_dt_(this->particles_->diffusion_dt_),
		  species_s_(species_s)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class FirstStageType>
	void SecondStageRK2<FirstStageType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < this->species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			all_species_[k][particle_i] = 0.5 * species_s_[m][particle_i] +
										  0.5 * (all_species_[k][particle_i] + dt * diffusion_dt_[m][particle_i]);
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	RelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::
		RelaxationOfAllDiffusionSpeciesRK2(typename FirstStageType::BodyRelationType &body_relation)
		: BaseDynamics<void>(body_relation.getSPHBody()),
		  rk2_initialization_(body_relation.getSPHBody(), species_s_),
		  rk2_1st_stage_(body_relation), rk2_2nd_stage_(body_relation, species_s_)
	{
		StdVec<BaseDiffusion *> species_diffusion_ =
			rk2_1st_stage_.diffusion_reaction_material_.SpeciesDiffusion();

		size_t number_of_diffusion_species = species_diffusion_.size();
		species_s_.resize(number_of_diffusion_species);

		constexpr int type_index = DataTypeIndex<Real>::value;
		for (size_t m = 0; m < number_of_diffusion_species; ++m)
		{
			// the size should be the same as that in the base particles
			species_s_[m].resize(rk2_1st_stage_.getParticles()->real_particles_bound_);
			// register data in base particles
			std::get<type_index>(rk2_1st_stage_.getParticles()->all_particle_data_).push_back(&species_s_[m]);
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	void RelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::exec(Real dt)
	{
		rk2_initialization_.exec();
		rk2_1st_stage_.exec(dt);
		rk2_2nd_stage_.exec(dt);
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	BaseRelaxationOfAllReactions<DiffusionReactionParticlesType>::
		BaseRelaxationOfAllReactions(SPHBody &sph_body)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<DiffusionReactionParticlesType>(sph_body),
		  reactive_species_(this->particles_->reactive_species_),
		  species_reaction_(this->particles_->diffusion_reaction_material_.SpeciesReaction()) {}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void BaseRelaxationOfAllReactions<DiffusionReactionParticlesType>::
		loadLocalSpecies(LocalSpecies &local_species, size_t index_i)
	{
		for (size_t k = 0; k != NumReactiveSpecies; ++k)
		{
			local_species[k] = (*reactive_species_[k])[index_i];
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void BaseRelaxationOfAllReactions<DiffusionReactionParticlesType>::
		applyGlobalSpecies(LocalSpecies &local_species, size_t index_i)
	{
		for (size_t k = 0; k != NumReactiveSpecies; ++k)
		{
			(*reactive_species_[k])[index_i] = local_species[k];
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void BaseRelaxationOfAllReactions<DiffusionReactionParticlesType>::
		advanceForwardStep(size_t index_i, Real dt)
	{
		LocalSpecies local_species;
		loadLocalSpecies(local_species, index_i);
		for (size_t k = 0; k != NumReactiveSpecies; ++k)
		{
			Real production_rate = species_reaction_.get_production_rates_[k](local_species);
			Real loss_rate = species_reaction_.get_loss_rates_[k](local_species);
			local_species[k] = updateAReactionSpecies(local_species[k], production_rate, loss_rate, dt);
		}
		applyGlobalSpecies(local_species, index_i);
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void BaseRelaxationOfAllReactions<DiffusionReactionParticlesType>::
		advanceBackwardStep(size_t index_i, Real dt)
	{
		LocalSpecies local_species;
		loadLocalSpecies(local_species, index_i);
		for (size_t k = NumReactiveSpecies; k != 0; --k)
		{
			size_t m = k - 1;
			Real production_rate = species_reaction_.get_production_rates_[m](local_species);
			Real loss_rate = species_reaction_.get_loss_rates_[m](local_species);
			local_species[m] = updateAReactionSpecies(local_species[m], production_rate, loss_rate, dt);
		}
		applyGlobalSpecies(local_species, index_i);
	}
	//=================================================================================================//
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP