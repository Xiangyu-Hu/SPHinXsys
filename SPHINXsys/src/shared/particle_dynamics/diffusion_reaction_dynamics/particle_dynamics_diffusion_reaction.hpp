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
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	DiffusionReactionInitialCondition<BodyType, BaseParticlesType, BaseMaterialType>::
		DiffusionReactionInitialCondition(BodyType &body)
		: ParticleDynamicsSimple(body),
		  DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(body),
		  pos_n_(this->particles_->pos_n_), species_n_(this->particles_->species_n_) {}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	GetDiffusionTimeStepSize<BodyType, BaseParticlesType, BaseMaterialType>::
		GetDiffusionTimeStepSize(BodyType &body)
		: ParticleDynamics<Real>(body),
		  DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(body)
	{
		Real smoothing_length = body.sph_adaptation_->ReferenceSmoothingLength();
		diff_time_step_ = this->material_->getDiffusionTimeStepSize(smoothing_length);
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::
		RelaxationOfAllDiffussionSpeciesInner(BaseBodyRelationInner &inner_relation)
		: InteractionDynamicsWithUpdate(*inner_relation.sph_body_),
		  DiffusionReactionInnerData<BodyType, BaseParticlesType, BaseMaterialType>(inner_relation),
		  species_n_(this->particles_->species_n_),
		  diffusion_dt_(this->particles_->diffusion_dt_), Vol_(this->particles_->Vol_)
	{
		species_diffusion_ = this->material_->SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::
		initializeDiffusionChangeRate(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			diffusion_dt_[m][particle_i] = 0;
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::
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
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_n_[k][particle_i] += dt * diffusion_dt_[m][particle_i];
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::
		Interaction(size_t index_i, Real dt)
	{
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType> *particles = this->particles_;
		Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];

		initializeDiffusionChangeRate(index_i);
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real dW_ij_ = inner_neighborhood.dW_ij_[n];
			Real r_ij_ = inner_neighborhood.r_ij_[n];
			Vecd &e_ij = inner_neighborhood.e_ij_[n];

			const Vecd &grad_ij = particles->getKernelGradient(index_i, index_j, dW_ij_, e_ij);
			Real area_ij = 2.0 * Vol_[index_j] * dot(grad_ij, e_ij) / r_ij_;
			getDiffusionChangeRate(index_i, index_j, e_ij, area_ij);
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::
		Update(size_t index_i, Real dt)
	{
		updateSpeciesDiffusion(index_i, dt);
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType,
			  class ContactBodyType, class ContactBaseParticlesType, class ContactBaseMaterialType>
	RelaxationOfAllDiffussionSpeciesComplex<BodyType, BaseParticlesType, BaseMaterialType,
											ContactBodyType, ContactBaseParticlesType, ContactBaseMaterialType>::
		RelaxationOfAllDiffussionSpeciesComplex(ComplexBodyRelation &complex_relation)
		: RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>(complex_relation.inner_relation_),
		  DiffusionReactionContactData<BodyType, BaseParticlesType, BaseMaterialType,
									   ContactBodyType, ContactBaseParticlesType, ContactBaseMaterialType>(complex_relation.contact_relation_),
		  species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_)
	{
		species_diffusion_ = this->material_->SpeciesDiffusion();

		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
			contact_species_n_.push_back(&(this->contact_particles_[k]->species_n_));
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType,
			  class ContactBodyType, class ContactBaseParticlesType, class ContactBaseMaterialType>
	void RelaxationOfAllDiffussionSpeciesComplex<BodyType, BaseParticlesType, BaseMaterialType,
												 ContactBodyType, ContactBaseParticlesType, ContactBaseMaterialType>::
		getDiffusionChangeRateContact(size_t particle_i, size_t particle_j, Vecd &e_ij,
									  Real surface_area_ij, StdVec<StdLargeVec<Real>> &species_n_k)
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
	template <class BodyType, class BaseParticlesType, class BaseMaterialType,
			  class ContactBodyType, class ContactBaseParticlesType, class ContactBaseMaterialType>
	void RelaxationOfAllDiffussionSpeciesComplex<BodyType, BaseParticlesType, BaseMaterialType,
												 ContactBodyType, ContactBaseParticlesType, ContactBaseMaterialType>::
		Interaction(size_t index_i, Real dt)
	{
		RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::Interaction(index_i, dt);
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
				Real dW_ij_ = contact_neighborhood.dW_ij_[n];
				Vecd &e_ij = contact_neighborhood.e_ij_[n];

				const Vecd &grad_ij = particles->getKernelGradient(index_i, index_j, dW_ij_, e_ij);
				Real area_ij = 2.0 * Vol_k[index_j] * dot(grad_ij, e_ij) / r_ij_;
				getDiffusionChangeRateContact(index_i, index_j, e_ij, area_ij, species_n_k);
			}
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	RungeKuttaInitialization<BodyType, BaseParticlesType, BaseMaterialType>::
		RungeKuttaInitialization(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &species_s)
		: ParticleDynamicsSimple(sph_body),
		  DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(sph_body),
		  species_n_(this->particles_->species_n_), species_s_(species_s)
	{
		species_diffusion_ = this->material_->SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RungeKuttaInitialization<BodyType, BaseParticlesType, BaseMaterialType>::
		initializeIntermediateValue(size_t particle_i)
	{
		for (size_t m = 0; m < species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_s_[m][particle_i] = species_n_[k][particle_i];
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RungeKuttaInitialization<BodyType, BaseParticlesType, BaseMaterialType>::
		Update(size_t index_i, Real dt)
	{
		initializeIntermediateValue(index_i);
	}
	//=================================================================================================//
	template <class RungeKutta2Stages1stStageType, class BodyRelationType>
	RungeKutta2Stages2ndStage<RungeKutta2Stages1stStageType, BodyRelationType>::
		RungeKutta2Stages2ndStage(BodyRelationType &body_relation, StdVec<StdLargeVec<Real>> &species_s)
		: RungeKutta2Stages1stStageType(body_relation),
		  species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_),
		  species_s_(species_s)
	{
		species_diffusion_ = this->material_->SpeciesDiffusion();
	}
	//=================================================================================================//
	template <class RungeKutta2Stages1stStageType, class BodyRelationType>
	void RungeKutta2Stages2ndStage<RungeKutta2Stages1stStageType, BodyRelationType>::
		updateSpeciesDiffusion(size_t particle_i, Real dt)
	{
		for (size_t m = 0; m < this->species_diffusion_.size(); ++m)
		{
			size_t k = species_diffusion_[m]->diffusion_species_index_;
			species_n_[k][particle_i] = 0.5 * species_s_[m][particle_i] + 0.5 * (species_n_[k][particle_i] + dt * diffusion_dt_[m][particle_i]);
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType,
			  class RungeKutta2Stages1stStageType, class BodyRelationType>
	RelaxationOfAllDiffusionSpeciesRK2<BodyType, BaseParticlesType, BaseMaterialType,
									   RungeKutta2Stages1stStageType, BodyRelationType>::
		RelaxationOfAllDiffusionSpeciesRK2(BodyRelationType &body_relation)
		: ParticleDynamics<void>(*body_relation.sph_body_),
		  DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(*body_relation.sph_body_),
		  runge_kutta_initialization_(*body_relation.sph_body_, species_s_),
		  runge_kutta_1st_stage_(body_relation),
		  runge_kutta_2nd_stage_(body_relation, species_s_)
	{
		StdVec<BaseDiffusion *> species_diffusion_ = this->material_->SpeciesDiffusion();

		size_t number_of_diffusion_species = species_diffusion_.size();
		species_s_.resize(number_of_diffusion_species);
		for (size_t m = 0; m < number_of_diffusion_species; ++m)
		{
			//the size should be the same as that in the base particles
			species_s_[m].resize(this->particles_->real_particles_bound_);
			//register data in base particles
			std::get<0>(this->particles_->all_particle_data_).push_back(&species_s_[m]);
		}
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType,
			  class RungeKutta2Stages1stStageType, class BodyRelationType>
	void RelaxationOfAllDiffusionSpeciesRK2<BodyType, BaseParticlesType, BaseMaterialType,
											RungeKutta2Stages1stStageType, BodyRelationType>::exec(Real dt)
	{
		runge_kutta_initialization_.exec();
		runge_kutta_1st_stage_.exec(dt);
		runge_kutta_2nd_stage_.exec(dt);
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType,
			  class RungeKutta2Stages1stStageType, class BodyRelationType>
	void RelaxationOfAllDiffusionSpeciesRK2<BodyType, BaseParticlesType, BaseMaterialType,
											RungeKutta2Stages1stStageType, BodyRelationType>::parallel_exec(Real dt)
	{
		runge_kutta_initialization_.parallel_exec();
		runge_kutta_1st_stage_.parallel_exec(dt);
		runge_kutta_2nd_stage_.parallel_exec(dt);
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	RelaxationOfAllReactionsForward<BodyType, BaseParticlesType, BaseMaterialType>::
		RelaxationOfAllReactionsForward(BodyType &body)
		: ParticleDynamicsSimple(body),
		  DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(body),
		  species_n_(this->particles_->species_n_)
	{
		species_reaction_ = this->material_->SpeciesReaction();
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllReactionsForward<BodyType, BaseParticlesType, BaseMaterialType>::
		Update(size_t index_i, Real dt)
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
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	RelaxationOfAllReactionsBackward<BodyType, BaseParticlesType, BaseMaterialType>::
		RelaxationOfAllReactionsBackward(BodyType &body)
		: ParticleDynamicsSimple(body),
		  DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(body),
		  species_n_(this->particles_->species_n_)
	{
		species_reaction_ = this->material_->SpeciesReaction();
	}
	//=================================================================================================//
	template <class BodyType, class BaseParticlesType, class BaseMaterialType>
	void RelaxationOfAllReactionsBackward<BodyType, BaseParticlesType, BaseMaterialType>::
		Update(size_t index_i, Real dt)
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
#endif //PARTICLE_DYNAMICS_DIFFUSION_REACTION_HPP