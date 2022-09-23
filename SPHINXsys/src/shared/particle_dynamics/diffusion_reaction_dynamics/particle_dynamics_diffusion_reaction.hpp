/**
 * @file 	particle_dynamics_diffusion_reaction.hpp
 * @author	Xiaojing Tang, Chi ZHang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP
#define PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP

#include "particle_dynamics_diffusion_reaction.h"

namespace SPH
{
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	DiffusionReactionInitialCondition<BaseParticlesType, BaseMaterialType>::
		DiffusionReactionInitialCondition(SPHBody &sph_body)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType>(sph_body),
		  pos_(this->particles_->pos_), species_n_(this->particles_->species_n_) {}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	GetDiffusionTimeStepSize<BaseParticlesType, BaseMaterialType>::
		GetDiffusionTimeStepSize(SPHBody &sph_body)
		: BaseDynamics<Real>(),
		  DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType>(sph_body)
	{
		Real smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
		diff_time_step_ = this->particles_->diffusion_reaction_material_.getDiffusionTimeStepSize(smoothing_length);
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>::
		RelaxationOfAllDiffusionSpeciesInner(BaseBodyRelationInner &inner_relation)
		: LocalDynamics(inner_relation.sph_body_),
		  DiffusionReactionInnerData<BaseParticlesType, BaseMaterialType>(inner_relation),
		  diffusion_reaction_material_(this->particles_->diffusion_reaction_material_),
		  species_n_(this->particles_->species_n_),
		  diffusion_dt_(this->particles_->diffusion_dt_), Vol_(this->particles_->Vol_)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>::
		initializeDiffusionChangeRate(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			diffusion_dt_[m][particle_i] = 0;
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>::
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
	template <class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_n_[k][particle_i] += dt * diffusion_dt_[m][particle_i];
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>::
		interaction(size_t index_i, Real dt)
	{
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType> *particles = this->particles_;
		Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];

		initializeDiffusionChangeRate(index_i);
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
			Real r_ij_ = inner_neighborhood.r_ij_[n];
			Vecd &e_ij = inner_neighborhood.e_ij_[n];

			const Vecd &grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
			Real area_ij = 2.0 * dot(grad_ijV_j, e_ij) / r_ij_;
			getDiffusionChangeRate(index_i, index_j, e_ij, area_ij);
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>::
		update(size_t index_i, Real dt)
	{
		updateSpeciesDiffusion(index_i, dt);
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
			  class ContactBaseParticlesType, class ContactBaseMaterialType>
	RelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
										   ContactBaseParticlesType, ContactBaseMaterialType>::
		RelaxationOfAllDiffusionSpeciesComplex(ComplexBodyRelation &complex_relation)
		: RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>(complex_relation.inner_relation_),
		  DiffusionReactionContactData<BaseParticlesType, BaseMaterialType,
									   ContactBaseParticlesType, ContactBaseMaterialType>(complex_relation.contact_relation_),
		  species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();

		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
			contact_species_n_.push_back(&(this->contact_particles_[k]->species_n_));
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
			  class ContactBaseParticlesType, class ContactBaseMaterialType>
	void RelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
												ContactBaseParticlesType, ContactBaseMaterialType>::
		getDiffusionChangeRateContact(size_t particle_i, size_t particle_j, Vecd &e_ij,
									  Real surface_area_ij, const StdVec<StdLargeVec<Real>> &species_n_k)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			Real diff_coff_ij = species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
			size_t l = species_diffusion_[m]->gradient_species_index_;
			Real phi_ij = species_n_[l][particle_i] - species_n_k[l][particle_j];
			diffusion_dt_[m][particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType,
			  class ContactBaseParticlesType, class ContactBaseMaterialType>
	void RelaxationOfAllDiffusionSpeciesComplex<BaseParticlesType, BaseMaterialType,
												ContactBaseParticlesType, ContactBaseMaterialType>::
		interaction(size_t index_i, Real dt)
	{
		RelaxationOfAllDiffusionSpeciesInner<BaseParticlesType, BaseMaterialType>::interaction(index_i, dt);
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType> *particles = this->particles_;

		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
			StdVec<StdLargeVec<Real>> &species_n_k = *(contact_species_n_[k]);

			Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
				size_t index_j = contact_neighborhood.j_[n];
				Real r_ij_ = contact_neighborhood.r_ij_[n];
				Real dW_ijV_j_ = contact_neighborhood.dW_ijV_j_[n];
				Vecd &e_ij = contact_neighborhood.e_ij_[n];

				const Vecd &grad_ijV_j = particles->getKernelGradient(index_i, index_j, dW_ijV_j_, e_ij);
				Real area_ij = 2.0 * dot(grad_ijV_j, e_ij) / r_ij_;
				getDiffusionChangeRateContact(index_i, index_j, e_ij, area_ij, species_n_k);
			}
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	InitializationRK<BaseParticlesType, BaseMaterialType>::
		InitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &species_s)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType>(sph_body),
		  species_n_(this->particles_->species_n_), species_s_(species_s)
	{
		species_diffusion_ = this->particles_->diffusion_reaction_material_.SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void InitializationRK<BaseParticlesType, BaseMaterialType>::
		initializeIntermediateValue(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_s_[m][particle_i] = species_n_[k][particle_i];
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void InitializationRK<BaseParticlesType, BaseMaterialType>::
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
		  species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_),
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
			species_n_[k][particle_i] = 0.5 * species_s_[m][particle_i] +
										0.5 * (species_n_[k][particle_i] + dt * diffusion_dt_[m][particle_i]);
		}
	}
	//=================================================================================================//
	template <class FirstStageType>
	RelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::
		RelaxationOfAllDiffusionSpeciesRK2(typename FirstStageType::BodyRelationType &body_relation)
		: BaseDynamics<void>(), rk2_initialization_(body_relation.sph_body_, species_s_),
		  rk2_1st_stage_(body_relation), rk2_2nd_stage_(body_relation, species_s_)
	{
		StdVec<BaseDiffusion *> species_diffusion_ = rk2_1st_stage_.diffusion_reaction_material_.SpeciesDiffusion();

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
	void RelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::exec(Real dt)
	{
		rk2_initialization_.exec();
		rk2_1st_stage_.exec(dt);
		rk2_2nd_stage_.exec(dt);
	}
	//=================================================================================================//
	template <class FirstStageType>
	void RelaxationOfAllDiffusionSpeciesRK2<FirstStageType>::parallel_exec(Real dt)
	{
		rk2_initialization_.parallel_exec();
		rk2_1st_stage_.parallel_exec(dt);
		rk2_2nd_stage_.parallel_exec(dt);
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	RelaxationOfAllReactionsForward<BaseParticlesType, BaseMaterialType>::
		RelaxationOfAllReactionsForward(SPHBody &sph_body)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType>(sph_body),
		  species_n_(this->particles_->species_n_)
	{
		species_reaction_ = this->particles_->diffusion_reaction_material_.SpeciesReaction();
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllReactionsForward<BaseParticlesType, BaseMaterialType>::
		update(size_t index_i, Real dt)
	{
		IndexVector &reactive_species = species_reaction_->reactive_species_;

		for (size_t m = 0; m != reactive_species.size(); ++m)
		{
			size_t k = reactive_species[m];
			Real production_rate = species_reaction_->get_production_rates_[k](species_n_, index_i);
			Real loss_rate = species_reaction_->get_loss_rates_[k](species_n_, index_i);
			species_n_[k][index_i] = updateAReactionSpecies(species_n_[k][index_i], production_rate, loss_rate, dt);
		}
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	RelaxationOfAllReactionsBackward<BaseParticlesType, BaseMaterialType>::
		RelaxationOfAllReactionsBackward(SPHBody &sph_body)
		: LocalDynamics(sph_body),
		  DiffusionReactionSimpleData<BaseParticlesType, BaseMaterialType>(sph_body),
		  species_n_(this->particles_->species_n_)
	{
		species_reaction_ = this->particles_->diffusion_reaction_material_.SpeciesReaction();
	}
	//=================================================================================================//
	template <class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllReactionsBackward<BaseParticlesType, BaseMaterialType>::
		update(size_t index_i, Real dt)
	{
		IndexVector &reactive_species = species_reaction_->reactive_species_;

		for (size_t m = reactive_species.size(); m != 0; --m)
		{
			size_t k = reactive_species[m - 1];
			Real production_rate = species_reaction_->get_production_rates_[k](species_n_, index_i);
			Real loss_rate = species_reaction_->get_loss_rates_[k](species_n_, index_i);
			species_n_[k][index_i] = updateAReactionSpecies(species_n_[k][index_i], production_rate, loss_rate, dt);
		}
	}
	//=================================================================================================//
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP