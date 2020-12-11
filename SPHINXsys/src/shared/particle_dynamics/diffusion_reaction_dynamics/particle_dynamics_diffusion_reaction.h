/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	particle_dynamics_diffusion.h
* @brief 	This is the particle dynamics aplliable for all type bodies
* @author	Chi ZHang and Xiangyu Hu
* @version	0.2
* 			Little knowledge of C++
*			You can do this : int a; const Real & b = a;
*			You get error when you do this : int a; Real & b = a; 
*			Note that this works fine for Visual Studio.
* 			a is of type int and is being converted to Real. 
*			So a temporary is created. Same is the case for user-defined types as well: Foo &obj = Foo(); 
*			See : https://stackoverflow.com/questions/18565167/non-const-lvalue-references
*			Chi Zhang
*/

#pragma once

#include "all_particle_dynamics.h"
#include "diffusion_reaction_particles.h"
#include "diffusion_reaction.h"


namespace SPH
{
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using  DiffusionReactionDataDelegateSimple = DataDelegateSimple<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>,
		DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionReactionDataDelegateInner = DataDelegateInner<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>,
		DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType,
		class ContactBodyType, class ContactBaseParticlesType, class ContactBaseMaterialType>
	using DiffusionReactionDataDelegateContact = DataDelegateContact<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>,
		DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>,
		ContactBodyType, DiffusionReactionParticles<ContactBaseParticlesType, ContactBaseMaterialType>,
		ContactBaseMaterialType, DataDelegateEmptyBase>;

	/**
	* @class  DiffusionReactionInitialCondition
	* @brief Computing the acoustic time step size
	* computing time step size
	 */
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class DiffusionReactionInitialCondition :
		public  ParticleDynamicsSimple,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	public:
		DiffusionReactionInitialCondition(BodyType* diffusion_body) :
			ParticleDynamicsSimple(diffusion_body),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(diffusion_body),
			pos_n_(this->particles_->pos_n_), species_n_(this->particles_->species_n_) {};
		virtual ~DiffusionReactionInitialCondition() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		StdVec<StdLargeVec<Real>>& species_n_;
	};

	/**
	* @class GetDiffusionTimeStepSize
	* @brief Computing the time step size based on diffusion coefficient and particle smoothing length
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class GetDiffusionTimeStepSize :
		public  ParticleDynamics<Real>,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	public:
		explicit GetDiffusionTimeStepSize(BodyType* body) :
			ParticleDynamics<Real>(body),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(body)
		{
			Real smoothing_length = body->kernel_->GetSmoothingLength();
			diff_time_step_ = this->material_->getDiffusionTimeStepSize(smoothing_length);
		};
		virtual ~GetDiffusionTimeStepSize() {};

		virtual Real exec(Real dt = 0.0) override { return diff_time_step_; };
		virtual Real parallel_exec(Real dt = 0.0) override { return exec(dt); };
	protected:
		Real diff_time_step_;
	};

	/**
	* @class RelaxationOfAllDiffussionSpeciesInner
	* @brief Compute the diffusion relaxation process of all species
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllDiffussionSpeciesInner :
		public InteractionDynamicsWithUpdate,
		public DiffusionReactionDataDelegateInner<BodyType, BaseParticlesType, BaseMaterialType>
	{
		/** all diffusion species and diffusion relation. */
		StdVec<BaseDiffusion*> species_diffusion_;
		StdVec<StdLargeVec<Real>>& species_n_;
		StdVec<StdLargeVec<Real>>& diffusion_dt_;
		StdLargeVec<Real>& Vol_;
	protected:
		/**
		  * @brief Initialize change rate to zero for all diffusion species.
		  * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i.
		  */
		void initializeDiffusionChangeRate(size_t particle_i)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				diffusion_dt_[m][particle_i] = 0;
			}
		};

		/**
		 * @brief Get change rate for all diffusion species.
		 * @param[in] particle_i Particle Index;
		 * @param[in] particle_j Particle Index;
		 * @param[in] e_ij Norm vector pointing from i to j;
		 * @param[in] surface_area_ij Surface area of particle interaction
		 * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i;
		 * @param[in] diffusion_reaction_data_j Diffusion reaction data of particle j;
		 */
		void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				Real diff_coff_ij = species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
				size_t l = species_diffusion_[m]->gradient_species_index_;
				Real phi_ij = species_n_[l][particle_i] - species_n_[l][particle_j];
				diffusion_dt_[m][particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
			}
		};

		/**
		 * @brief Update all diffusion species.
		 * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i;
		 * @param[in] dt Time step;
		 */
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				species_n_[k][particle_i] += dt * diffusion_dt_[m][particle_i];
			}
		};

		virtual void Interaction(size_t index_i, Real dt = 0.0) override
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];

			initializeDiffusionChangeRate(index_i);
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij_ = inner_neighborhood.dW_ij_[n];
				Real r_ij_ = inner_neighborhood.r_ij_[n];
				Vecd& e_ij = inner_neighborhood.e_ij_[n];

				const Vecd& gradi_ij = particles->getKernelGradient(index_i, index_j, dW_ij_, e_ij);
				Real area_ij = 2.0 * Vol_[index_j] * dot(gradi_ij, e_ij) / r_ij_;
				getDiffusionChangeRate(index_i, index_j, e_ij, area_ij);
			}
		};

		virtual void Update(size_t index_i, Real dt = 0.0) override
		{
			updateSpeciesDiffusion(index_i, dt);
		};
	public:
		RelaxationOfAllDiffussionSpeciesInner(SPHBodyInnerRelation* body_inner_relation)
			: InteractionDynamicsWithUpdate(body_inner_relation->sph_body_),
			DiffusionReactionDataDelegateInner<BodyType, BaseParticlesType, BaseMaterialType>(body_inner_relation),
			species_n_(this->particles_->species_n_),
			diffusion_dt_(this->particles_->diffusion_dt_), Vol_(this->particles_->Vol_)
		{
			species_diffusion_ = this->material_->SpeciesDiffusion();
		};
		virtual ~RelaxationOfAllDiffussionSpeciesInner() {};
	};

	/**
	 * @class RelaxationOfAllDiffussionSpeciesComplex
	 *Complex diffusion relaxation between two different bodies
	 */
	template<class BodyType, class BaseParticlesType, class BaseMaterialType,
		class ContactBodyType, class ContactBaseParticlesType, class ContactBaseMaterialType>
		class RelaxationOfAllDiffussionSpeciesComplex :
		public RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>,
		public DiffusionReactionDataDelegateContact<BodyType, BaseParticlesType, BaseMaterialType,
		ContactBodyType, ContactBaseParticlesType, ContactBaseMaterialType>
	{
		/** all diffusion species and diffusion relation. */
		StdVec<BaseDiffusion*> species_diffusion_;
		StdVec<StdLargeVec<Real>>& species_n_;
		StdVec<StdLargeVec<Real>>& diffusion_dt_;
	protected:
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec<StdVec< StdLargeVec<Real>>*> contact_species_n_;

		//get diffusion change rate for the species  contact  
		void getDiffusionChangeRateContact(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij, StdVec<StdLargeVec<Real>>& species_n_k)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				Real diff_coff_ij = species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
				size_t l = species_diffusion_[m]->gradient_species_index_;
				Real phi_ij = species_n_[l][particle_i] - species_n_k[l][particle_j];
				diffusion_dt_[m][particle_i] += diff_coff_ij * phi_ij * surface_area_ij;
			}
		};

		//complex interaction between two contact bodies
		virtual void Interaction(size_t index_i, Real dt = 0.0) override
		{
			//inner interaction
			RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>::Interaction(index_i, dt);
			/** Contact interaction. */
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				StdVec<StdLargeVec<Real>>& species_n_k = *(contact_species_n_[k]);

				Neighborhood& contact_neighborhood = this->contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Real r_ij_ = contact_neighborhood.r_ij_[n];
					Real dW_ij_ = contact_neighborhood.dW_ij_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];

					const Vecd& gradi_ij = particles->getKernelGradient(index_i, index_j, dW_ij_, e_ij);
					Real area_ij = 2.0 * Vol_k[index_j] * dot(gradi_ij, e_ij) / r_ij_;
					getDiffusionChangeRateContact(index_i, index_j, e_ij, area_ij, species_n_k);
				}
			}
		};

	public:
		RelaxationOfAllDiffussionSpeciesComplex(SPHBodyComplexRelation* body_complex_relation) :
			RelaxationOfAllDiffussionSpeciesInner<BodyType, BaseParticlesType, BaseMaterialType>(body_complex_relation->InnerRelation()),
			DiffusionReactionDataDelegateContact<BodyType, BaseParticlesType, BaseMaterialType,
			ContactBodyType, ContactBaseParticlesType, ContactBaseMaterialType>(body_complex_relation->ContactRelation()),
			species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_)
		{
			species_diffusion_ = this->material_->SpeciesDiffusion();

			for (size_t k = 0; k != this->contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
				contact_species_n_.push_back(&(this->contact_particles_[k]->species_n_));
			}
		};
		virtual ~RelaxationOfAllDiffussionSpeciesComplex() {};

	};

	/**
	* @class RungeKuttaInitialization
	* @brief initialization of a runge-kutta integration scheme
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RungeKuttaInitialization :
		public ParticleDynamicsSimple,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
		StdVec<BaseDiffusion*> species_diffusion_;
		StdVec<StdLargeVec<Real>>& species_n_, & species_s_;

		void initializeIntermediateValue(size_t particle_i)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				species_s_[m][particle_i] = species_n_[k][particle_i];
			}
		};

		virtual void Update(size_t index_i, Real dt = 0.0) override
		{
			initializeIntermediateValue(index_i);
		};
	public:
		RungeKuttaInitialization(SPHBody* body, StdVec<StdLargeVec<Real>>& species_s) :
			ParticleDynamicsSimple(body),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(body),
			species_n_(this->particles_->species_n_), species_s_(species_s)
		{
			species_diffusion_ = this->material_->SpeciesDiffusion();
		};
		virtual ~RungeKuttaInitialization() {};
	};

	/**
	* @class RungeKutta2Stages2ndStage
	* @brief the second stage of the second runge-kutta scheme
	*/
	template<class RungeKutta2Stages1stStageType, class BodyRelationType>
	class RungeKutta2Stages2ndStage : public RungeKutta2Stages1stStageType
	{
		/** all diffusion species and diffusion relation. */
		StdVec<BaseDiffusion*> species_diffusion_;
		StdVec<StdLargeVec<Real>>& species_n_;
		StdVec<StdLargeVec<Real>>& diffusion_dt_;
	protected:
		StdVec<StdLargeVec<Real>>& species_s_;
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override
		{
			for (size_t m = 0; m < this->species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				species_n_[k][particle_i] = 0.5 * species_s_[m][particle_i]
					+ 0.5 * (species_n_[k][particle_i] + dt * diffusion_dt_[m][particle_i]);
			}
		};
	public:
		RungeKutta2Stages2ndStage(BodyRelationType* body_relation, StdVec<StdLargeVec<Real>>& species_s)
			: RungeKutta2Stages1stStageType(body_relation),
			species_n_(this->particles_->species_n_), diffusion_dt_(this->particles_->diffusion_dt_),
			species_s_(species_s)
		{
			species_diffusion_ = this->material_->SpeciesDiffusion();
		};
		virtual ~RungeKutta2Stages2ndStage() {};
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRK2
	 * @brief Compute the diffusion relaxation process of all species
	 * with second order Runge-Kutta time stepping
	 */
	template<class BodyType, class BaseParticlesType, class BaseMaterialType,
		class RungeKutta2Stages1stStageType, class BodyRelationType>
	class RelaxationOfAllDiffusionSpeciesRK2 : public ParticleDynamics<void>,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	protected:
		StdVec<BaseDiffusion*> species_diffusion_;
		/** Intermediate Value */
		StdVec<StdLargeVec<Real>> species_s_;

		RungeKuttaInitialization<BodyType, BaseParticlesType, BaseMaterialType> runge_kutta_initialization_;
		RungeKutta2Stages1stStageType runge_kutta_1st_stage_;
		RungeKutta2Stages2ndStage<RungeKutta2Stages1stStageType, BodyRelationType> runge_kutta_2nd_stage_;
	public:
		RelaxationOfAllDiffusionSpeciesRK2(BodyRelationType* body_relation) :
			ParticleDynamics<void>(body_relation->sph_body_),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(body_relation->sph_body_),
			runge_kutta_initialization_(body_relation->sph_body_, species_s_),
			runge_kutta_1st_stage_(body_relation),
			runge_kutta_2nd_stage_(body_relation, species_s_)
		{
			StdVec<BaseDiffusion*> species_diffusion_ = this->material_->SpeciesDiffusion();

			size_t number_of_diffusion_species = species_diffusion_.size();
			species_s_.resize(number_of_diffusion_species);
			for (size_t m = 0; m < number_of_diffusion_species; ++m)
			{
				//the size should be the same as that in the base particles
				species_s_[m].resize(this->particles_->real_particles_bound_);
				//register data in base particles
				this->particles_->registered_scalars_.push_back(&species_s_[m]);
			}
		};
		virtual ~RelaxationOfAllDiffusionSpeciesRK2() {};

		virtual void exec(Real dt = 0.0) override {
			runge_kutta_initialization_.exec();
			runge_kutta_1st_stage_.exec(dt);
			runge_kutta_2nd_stage_.exec(dt);
		};
		virtual void parallel_exec(Real dt = 0.0) override {
			runge_kutta_initialization_.parallel_exec();
			runge_kutta_1st_stage_.parallel_exec(dt);
			runge_kutta_2nd_stage_.parallel_exec(dt);
		};
	};

	/**
	* @class RelaxationOfAllReactionsForward
	* @brief Compute the reaction process of all species by forward splitting
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllReactionsForward :
		public ParticleDynamicsSimple,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
		/** The reaction model for all reactive species. */
		BaseReactionModel* species_reaction_;
		StdVec<StdLargeVec<Real>>& species_n_;
	protected:
		/**
		 * @brief Splitting scheme for directly computing one time step integration for a species.
		 * @param[in] input Input state of species.
		 * @param[in] production_rate Production rate of species.
		 * @param[in] loss_rate Loss rate of species.
		 * @param[in] dt Time step size
		 * @param[out] Change rate of species.
		 **/
		Real updateAReactionSpecies(Real input, Real production_rate, Real loss_rate, Real dt)
		{
			return input * exp(-loss_rate * dt) + production_rate * (1.0 - exp(-loss_rate * dt)) / (loss_rate + TinyReal);
		};
		/** Get change rate for all reactive species by forward sweeping.
		 * @brief Get change rate for all reactive species by backward sweeping.
		 * @param[in] diffusion_reaction_data_i Diffusion Reaction Data.
		 * @param[in] dt Time step size.
		 **/
		void UpdateReactiveSpeciesForward(size_t particle_i, Real dt) {

		};

		virtual void Update(size_t index_i, Real dt = 0.0) override
		{
			IndexVector& reactive_species = species_reaction_->reactive_species_;

			for (size_t m = 0; m != reactive_species.size(); ++m) {
				size_t k = reactive_species[m];
				Real production_rate = species_reaction_->get_production_rates_[k](species_n_, index_i);
				Real loss_rate = species_reaction_->get_loss_rates_[k](species_n_, index_i);
				Real input = species_n_[k][index_i];
				species_n_[k][index_i] = updateAReactionSpecies(input, production_rate, loss_rate, dt);
			}
		}
	public:
		RelaxationOfAllReactionsForward(BodyType* body) :
			ParticleDynamicsSimple(body),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(body),
			species_n_(this->particles_->species_n_)
		{
			species_reaction_ = this->material_->SpeciesReaction();
		};
		virtual ~RelaxationOfAllReactionsForward() {};
	};

	/**
	* @class RelaxationOfAllReactionsBackward
	* @brief Compute the reaction process of all species by backward splitting
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllReactionsBackward :
		public ParticleDynamicsSimple,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
		/** The reaction model for all reactive species. */
		BaseReactionModel* species_reaction_;
		StdVec<StdLargeVec<Real>>& species_n_;
	protected:
		/**
		 * @brief Splitting scheme for directly computing one time step integeration for a species.
		 * @param[in] input Input state of species.
		 * @param[in] production_rate Production rate of species.
		 * @param[in] loss_rate Loss rate of species.
		 * @param[in] dt Time step size
		 * @param[out] Change rate of species.
		 **/
		Real updateAReactionSpecies(Real input, Real production_rate, Real loss_rate, Real dt)
		{
			return input * exp(-loss_rate * dt) + production_rate * (1.0 - exp(-loss_rate * dt)) / (loss_rate + TinyReal);
		};

		virtual void Update(size_t index_i, Real dt = 0.0) override
		{
			IndexVector& reactive_species = species_reaction_->reactive_species_;

			for (size_t m = reactive_species.size(); m != 0; --m) {
				size_t k = reactive_species[m - 1];
				Real production_rate = species_reaction_->get_production_rates_[k](species_n_, index_i);
				Real loss_rate = species_reaction_->get_loss_rates_[k](species_n_, index_i);
				Real input = species_n_[k][index_i];
				species_n_[k][index_i] = updateAReactionSpecies(input, production_rate, loss_rate, dt);
			}
		};
	public:
		RelaxationOfAllReactionsBackward(BodyType* body) :
			ParticleDynamicsSimple(body),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(body),
			species_n_(this->particles_->species_n_)
		{
			species_reaction_ = this->material_->SpeciesReaction();
		};
		virtual ~RelaxationOfAllReactionsBackward() {};
	};

	/**
	 * @class DiffusionBoundaryCondition
	 * @brief set boundary condition for diffusion problem
	 */
	template <class BodyType, class BaseParticlesType, class BodyPartByParticleType, class BaseMaterialType>
	class DiffusionBoundaryCondition :
		public PartDynamicsByParticle,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	public:
		DiffusionBoundaryCondition(BodyType* body, BodyPartByParticleType* body_part) :
			PartDynamicsByParticle(body, body_part),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(body),
			pos_n_(this->particles_->pos_n_), species_n_(this->particles_->species_n_) {};
		virtual ~DiffusionBoundaryCondition() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		StdVec<StdLargeVec<Real>>& species_n_;

	};

	/**
	* @class DiffusionBasedMapping
	* @brief Mapping inside of body according to diffusion.
	* This is a abstract class to be override for case specific implementation
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class DiffusionBasedMapping :
		public ParticleDynamicsSimple,
		public DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	public:
		DiffusionBasedMapping(BodyType* body) :
			ParticleDynamicsSimple(body),
			DiffusionReactionDataDelegateSimple<BodyType, BaseParticlesType, BaseMaterialType>(body),
			pos_n_(this->particles_->pos_n_), species_n_(this->particles_->species_n_) {};
		virtual ~DiffusionBasedMapping() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		StdVec<StdLargeVec<Real>>& species_n_;
	};
}
