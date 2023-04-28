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
 *  HU1527/12-1 and HU1527/12-4												*
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
		  material_(this->particles_->diffusion_reaction_material_),
		  all_diffusions_(material_.AllDiffusions()),
		  diffusion_species_(this->particles_->DiffusionSpecies()),
		  gradient_species_(this->particles_->GradientSpecies())
	{
		diffusion_dt_.resize(all_diffusions_.size());
		StdVec<std::string> &all_species_names =this->particles_->AllSpeciesNames();
		IndexVector &diffusion_species_indexes = material_.DiffusionSpeciesIndexes();
		for (size_t i = 0; i != all_diffusions_.size(); ++i)
		{
			// Register specie change rate as shared variable
			std::string &diffusion_species_name = all_species_names[diffusion_species_indexes[i]];
			diffusion_dt_[i] = this->particles_->template registerSharedVariable<Real>(diffusion_species_name + "ChangeRate");
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		initializeDiffusionChangeRate(size_t particle_i)
	{
		for (size_t m = 0; m < all_diffusions_.size(); ++m)
		{
			(*diffusion_dt_[m])[particle_i] = 0;
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij)
	{
		for (size_t m = 0; m < all_diffusions_.size(); ++m)
		{
			Real diff_coff_ij =
				all_diffusions_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
			StdLargeVec<Real> &gradient_species = *gradient_species_[m];
			Real phi_ij = gradient_species[particle_i] - gradient_species[particle_j];
			(*diffusion_dt_[m])[particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < all_diffusions_.size(); ++m)
		{
			(*diffusion_species_[m])[particle_i] += dt * (*diffusion_dt_[m])[particle_i];
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
			  ContactDiffusionReactionParticlesType>(complex_relation.getContactRelation())
	{
		contact_gradient_species_.resize(this->contact_particles_.size());

		StdVec<std::string> &all_species_names = this->particles_->AllSpeciesNames();
		for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
		{
			size_t l = this->all_diffusions_[m]->gradient_species_index_;
			std::string &inner_species_name_m = all_species_names[l];
			for (size_t k = 0; k != this->contact_particles_.size(); ++k)
			{
				size_t contact_species_index_k_m = this->contact_particles_[k]->AllSpeciesIndexMap()[inner_species_name_m];
				StdVec<std::string> &all_contact_species_names_k = this->contact_particles_[k]->AllSpeciesNames();
				std::string &contact_species_name_k_m = all_contact_species_names_k[contact_species_index_k_m];

				if (inner_species_name_m != contact_species_name_k_m)
				{
					std::cout << "\n Error: inner species '" << inner_species_name_m
							  << "' and contact species '" << contact_species_name_k_m << "' not match! " << std::endl;
					std::cout << __FILE__ << ':' << __LINE__ << std::endl;
					exit(1);
				}
				StdVec<StdLargeVec<Real>> &all_contact_species_k = this->contact_particles_[k]->all_species_;

				contact_gradient_species_[k].push_back(&all_contact_species_k[contact_species_index_k_m]);
			}
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	void RelaxationOfAllDiffusionSpeciesComplex<
		DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>::
		getDiffusionChangeRateContact(size_t particle_i, size_t particle_j, Vecd &e_ij,
									  Real surface_area_ij, const StdVec<StdLargeVec<Real> *> &gradient_species_k)
	{
		for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
		{
			Real diff_coff_ij =
				this->all_diffusions_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
			Real phi_ij = (*this->diffusion_species_[m])[particle_i] - (*gradient_species_k[m])[particle_j];
			(*this->diffusion_dt_[m])[particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
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
			StdVec<StdLargeVec<Real> *> &gradient_species_k = contact_gradient_species_[k];

			Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real r_ij_ = contact_neighborhood.r_ij_[n];
				Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd &e_ij = contact_neighborhood.e_ij_[n];

				const Vecd &grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
				Real area_ij = 2.0 * grad_ijV_j.dot(e_ij) / r_ij_;
				getDiffusionChangeRateContact(index_i, index_j, e_ij, area_ij, gradient_species_k);
			}
		}
	}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	InitializationRK<DiffusionReactionParticlesType>::
		InitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &diffusion_species_s)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<DiffusionReactionParticlesType>(sph_body),
		  material_(this->particles_->diffusion_reaction_material_),
		  all_diffusions_(material_.AllDiffusions()),
		  diffusion_species_(this->particles_->DiffusionSpecies()),
		  diffusion_species_s_(diffusion_species_s) {}
	//=================================================================================================//
	template <class DiffusionReactionParticlesType>
	void InitializationRK<DiffusionReactionParticlesType>::
		update(size_t index_i, Real dt)
	{
		for (size_t m = 0; m < all_diffusions_.size(); ++m)
		{
			diffusion_species_s_[m][index_i] = (*diffusion_species_[m])[index_i];
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	SecondStageRK2<FirstStageType>::
		SecondStageRK2(typename FirstStageType::BodyRelationType &body_relation,
					   StdVec<StdLargeVec<Real>> &diffusion_species_s)
		: FirstStageType(body_relation), diffusion_species_s_(diffusion_species_s) {}
	//=================================================================================================//
	template <class FirstStageType>
	void SecondStageRK2<FirstStageType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
		{
			(*this->diffusion_species_[m])[particle_i] =
				0.5 * diffusion_species_s_[m][particle_i] +
				0.5 * ((*this->diffusion_species_[m])[particle_i] + dt * (*this->diffusion_dt_[m])[particle_i]);
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	RelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::
		RelaxationOfAllDiffusionSpeciesRK2(typename FirstStageType::BodyRelationType &body_relation)
		: BaseDynamics<void>(body_relation.getSPHBody()),
		  rk2_initialization_(body_relation.getSPHBody(), diffusion_species_s_),
		  rk2_1st_stage_(body_relation), rk2_2nd_stage_(body_relation, diffusion_species_s_),
		  all_diffusions_(rk2_1st_stage_.AllDiffusions())
	{
		diffusion_species_s_.resize(all_diffusions_.size());
		StdVec<std::string> &all_species_names = rk2_1st_stage_.getParticles()->AllSpeciesNames();
		for (size_t i = 0; i != all_diffusions_.size(); ++i)
		{
			// register diffusion species intermediate
			size_t diffusion_species_index = all_diffusions_[i]->diffusion_species_index_;
			std::string &diffusion_species_name = all_species_names[diffusion_species_index];
			rk2_1st_stage_.getParticles()->registerVariable(diffusion_species_s_[i], diffusion_species_name + "Intermediate");
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
		  reactive_species_(this->particles_->ReactiveSpecies()),
		  reaction_model_(this->particles_->diffusion_reaction_material_.ReactionModel()) {}
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
			Real production_rate = reaction_model_.get_production_rates_[k](local_species);
			Real loss_rate = reaction_model_.get_loss_rates_[k](local_species);
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
			Real production_rate = reaction_model_.get_production_rates_[m](local_species);
			Real loss_rate = reaction_model_.get_loss_rates_[m](local_species);
			local_species[m] = updateAReactionSpecies(local_species[m], production_rate, loss_rate, dt);
		}
		applyGlobalSpecies(local_species, index_i);
	}
	//=================================================================================================//
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP