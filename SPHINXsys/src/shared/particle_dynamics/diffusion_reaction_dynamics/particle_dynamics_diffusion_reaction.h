/**
* @file 	particle_dynamics_diffusion.h
* @brief 	This is the particle dynnamics apllicable for all type bodies
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
*			0.2.1
*			Low-storage 4rd Runge-Kutta algorithm is added.
*			Ref : doi:10.1016/j.jcp.2009.11.006
*			A template-parameter shall not be redeclared within its scope (including nested scopes). 
*			A template- parameter shall not have the same name as the template name. 
*/

#pragma once

#include "all_particle_dynamics.h"
#include "diffusion_reaction_particles.h"
#include "diffusion_reaction.h"


namespace SPH
{
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionBase = ParticleDynamics<Real, BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionReactionSimple = ParticleDynamicsSimple<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionInner = ParticleDynamicsInner<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionInnerWithUpdate = ParticleDynamicsInnerWithUpdate<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>>;

	template <class BodyType, class BaseParticlesType, class BodyPartByParticleType, class BaseMaterialType>
	using DiffusionReactionConstraint = ConstraintByParticle<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, BodyPartByParticleType, 
		DiffusionReactionMaterial<BaseParticlesType, BaseMaterialType>>;

	/**
	* @class GetDiffusionTimeStepSize
	* @brief Computing the acoustic time step size
	* computing time step size
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class GetDiffusionTimeStepSize
		: public  DiffusionBase<BodyType, BaseParticlesType, BaseMaterialType>
	{
	protected:
		Real diff_time_step_;
	public:
		explicit GetDiffusionTimeStepSize(BodyType* body)
			: DiffusionBase<BodyType, BaseParticlesType, BaseMaterialType>(body) 
		{
			Real smoothing_length = body->kernel_->GetSmoothingLength();
			diff_time_step_ = this->material_->getDiffusionTimeStepSize(smoothing_length);
		};
		virtual ~GetDiffusionTimeStepSize() {};

		virtual Real exec(Real dt = 0.0) override { return diff_time_step_; };
		virtual Real parallel_exec(Real dt = 0.0) override { return exec(dt); };
	};

	/**
	* @class RelaxationOfAllDifussionSpecies
	* @brief Compute the diffusion relaxation process of all species
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllDifussionSpecies
		: public DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>
	{
		/** all diffusion species and diffusion relation. */
		StdVec<BaseDiffusion*> species_diffusion_;
	protected:
		/**
		  * @brief Initialize change rate to zero for all diffusion species.
		  * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i.
		  */
		void initializeDiffusionChangeRate(DiffusionReactionData& diffusion_reaction_data_i)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				diffusion_reaction_data_i.dspecies_dt_[k] = 0;
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
		void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij,
			DiffusionReactionData& diffusion_reaction_data_i, DiffusionReactionData& diffusion_reaction_data_j)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				Real diff_coff_ij = species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				size_t l = species_diffusion_[m]->gradient_species_index_;
				Real phi_ij = diffusion_reaction_data_i.species_n_[k] - diffusion_reaction_data_j.species_n_[k];
				diffusion_reaction_data_i.dspecies_dt_[k] += diff_coff_ij * phi_ij * surface_area_ij;
			}
		};

		/**
		 * @brief Update all diffusion species.
		 * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i;
		 * @param[in] dt Time step;
		 */
		virtual void updateDiffusionSpecies(DiffusionReactionData& diffusion_reaction_data_i, Real dt)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				diffusion_reaction_data_i.species_n_[k] += dt * diffusion_reaction_data_i.dspecies_dt_[k];
			}
		};

		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			Neighborhood& neighborhood = (*this->inner_configuration_)[index_particle_i];

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			StdLargeVec<BaseParticleData>& base_particle_data = particles->base_particle_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];

			initializeDiffusionChangeRate(diffusion_reaction_data_i);
			NeighborList& neighors = std::get<0>(neighborhood);
			for (size_t n = 0; n != std::get<2>(neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real Vol_j = base_particle_data[index_particle_j].Vol_;
				DiffusionReactionData& diffusion_reaction_data_j = diffusion_reaction_data[index_particle_j];
	
				const Vecd& gradi_ij = particles->getKernelGradient(index_particle_i, index_particle_j, neighboring_particle->dW_ij_, e_ij);
				Real area_ij = 2.0 * base_particle_data[index_particle_j].Vol_ * dot(gradi_ij, e_ij) / neighboring_particle->r_ij_;
				getDiffusionChangeRate(index_particle_i, index_particle_j, e_ij, area_ij,
					diffusion_reaction_data_i, diffusion_reaction_data_j);
			}
		};

		virtual void Update(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			updateDiffusionSpecies(diffusion_reaction_data_i, dt);
		};
	public:
		RelaxationOfAllDifussionSpecies(BodyType* body)
			: DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>(body) {
			species_diffusion_ = this->material_->getDiffusionSpecies();
		};
		virtual ~RelaxationOfAllDifussionSpecies() {};
	};

	/**
		* @class RelaxationOfAllDifussionSpeciesRK2
		* @brief Compute the diffusion relaxation process of all species
		* with second order Runge Kutta time stepping
		*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllDifussionSpeciesRK2
		: public ParticleDynamics<void, BodyType, BaseParticlesType, BaseMaterialType>
	{
	protected:
		template<class ParameterBody, class ParameterBaseParticles, class ParameterBaseMaterial>
		class RungeKuttaInitialization :
			public DiffusionReactionSimple<ParameterBody, ParameterBaseParticles, ParameterBaseMaterial>
		{
			StdVec<BaseDiffusion*> species_diffusion_;
		protected:

			void initializeIntermediateValue(DiffusionReactionData& diffusion_reaction_data_i)
			{
				for (size_t m = 0; m < species_diffusion_.size(); ++m)
				{
					size_t k = species_diffusion_[m]->diffusion_species_index_;
					diffusion_reaction_data_i.species_s_[k] = diffusion_reaction_data_i.species_n_[k];
				}
			};

			virtual void Update(size_t index_particle_i, Real dt = 0.0) override 
			{
				DiffusionReactionParticles<ParameterBaseParticles, ParameterBaseMaterial>* particles = this->particles_;
				StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
				DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
				initializeIntermediateValue(diffusion_reaction_data_i);
			};
		public:
			RungeKuttaInitialization(ParameterBody* body)
				: DiffusionReactionSimple<ParameterBody, ParameterBaseParticles, ParameterBaseMaterial>(body) 
			{
				species_diffusion_ = this->material_->getDiffusionSpecies();
			};;
			virtual ~RungeKuttaInitialization() {};
		};

		template<class ParameterBody, class ParameterBaseParticles, class ParameterBaseMaterial>
		class RungeKutta2Stages2ndStage : public RelaxationOfAllDifussionSpecies<ParameterBody, ParameterBaseParticles, ParameterBaseMaterial>
		{
			StdVec<BaseDiffusion*> species_diffusion_;
		protected:
			virtual void updateDiffusionSpecies(DiffusionReactionData& diffusion_reaction_data_i, Real dt) override
			{
				for (size_t m = 0; m < species_diffusion_.size(); ++m)
				{
					size_t k = species_diffusion_[m]->diffusion_species_index_;
					diffusion_reaction_data_i.species_n_[k] = 0.5 * diffusion_reaction_data_i.species_s_[k]
					+ 0.5 * (diffusion_reaction_data_i.species_n_[k] + dt * diffusion_reaction_data_i.dspecies_dt_[k]);
				}
			};
		public:
			RungeKutta2Stages2ndStage(ParameterBody* body)
				: RelaxationOfAllDifussionSpecies<ParameterBody, ParameterBaseParticles, ParameterBaseMaterial>(body) 
			{
				species_diffusion_ = this->material_->getDiffusionSpecies();
			};
			virtual ~RungeKutta2Stages2ndStage() {};
		};

		RungeKuttaInitialization<BodyType, BaseParticlesType, BaseMaterialType> runge_kutta_initialization_;
		RelaxationOfAllDifussionSpecies<BodyType, BaseParticlesType, BaseMaterialType> runge_kutta_1st_stage_;
		RungeKutta2Stages2ndStage<BodyType, BaseParticlesType, BaseMaterialType> runge_kutta_2nd_stage_;
	public:
		RelaxationOfAllDifussionSpeciesRK2(BodyType* body)
			: ParticleDynamics<void, BodyType, BaseParticlesType, BaseMaterialType>(body),
			runge_kutta_initialization_(body), runge_kutta_1st_stage_(body),
			runge_kutta_2nd_stage_(body){};
		virtual ~RelaxationOfAllDifussionSpeciesRK2() {};

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
	* @class RelaxationOfAllDifussionSpeciesRungeKutta
	* @brief Compute the diffusion relaxation process of all species using 4th Runge-Kuttas cheme
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllDifussionSpeciesRungeKutta
		: public DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>
	{
	private:
		Real gamma_1_[5] = { 0.0, 0.0, 0.12109848,-3.8438337, 0.5463709 };
		Real gamma_2_[5] = { 0.0, 1.0, 0.7217817,  2.1212093, 0.1986530 };
		Real beta_[5] = { 0.0, 1.1937439, 0.0992799, 1.1316780, 0.3106657 };
		Real delta_[4] = { 1.0, 0.2176833, 1.0658413, 0.0 };
		int RK_step_;
	protected:
		/** all diffusion species and diffusion relation. */
		StdVec<BaseDiffusion*> species_diffusion_;
		/**
		 * @brief Initialize the stages for Runge-Kutta scheme.
		 * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i.
		 */
		void initializeStageForRungeKutta(DiffusionReactionData& diffusion_reaction_data_i)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				diffusion_reaction_data_i.species_s_[k] = 0.0;
			}
		};

		/**
		 * @brief Update all diffusion species.
		 * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i;
		 * @param[in] dt Time step;
		 */
		void updateStageforRungeKutta(DiffusionReactionData& diffusion_reaction_data_i, Real delta)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				diffusion_reaction_data_i.species_s_[k] += delta * diffusion_reaction_data_i.species_n_[k];
				diffusion_reaction_data_i.dspecies_dt_[k] = 0;
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
		void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij,
			DiffusionReactionData& diffusion_reaction_data_i, DiffusionReactionData& diffusion_reaction_data_j)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				Real diff_coff_ij = species_diffusion_[m]->getInterParticleDiffusionCoff(particle_i, particle_j, e_ij);
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				size_t l = species_diffusion_[m]->gradient_species_index_;
				Real phi_ij = diffusion_reaction_data_i.species_n_[k] - diffusion_reaction_data_j.species_n_[k];
				diffusion_reaction_data_i.dspecies_dt_[k] += diff_coff_ij * phi_ij * surface_area_ij;
			}
		};

		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			Neighborhood& neighborhood = (*this->inner_configuration_)[index_particle_i];

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			StdLargeVec<BaseParticleData>& base_particle_data = particles->base_particle_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			if (RK_step_ == 1)
			{
				initializeStageForRungeKutta(diffusion_reaction_data_i);
			}
			
			updateStageforRungeKutta(diffusion_reaction_data_i, delta_[RK_step_ - 1]);
			NeighborList& neighors = std::get<0>(neighborhood);
			for (size_t n = 0; n != std::get<2>(neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real Vol_j = base_particle_data[index_particle_j].Vol_;
				DiffusionReactionData& diffusion_reaction_data_j = diffusion_reaction_data[index_particle_j];

				const Vecd& gradi_ij = particles->getKernelGradient(index_particle_i, index_particle_j, neighboring_particle->dW_ij_, e_ij);
				Real area_ij = 2.0 * base_particle_data[index_particle_j].Vol_ * dot(gradi_ij, e_ij) / neighboring_particle->r_ij_;
				getDiffusionChangeRate(index_particle_i, index_particle_j, e_ij, area_ij,
					diffusion_reaction_data_i, diffusion_reaction_data_j);
			}
		};
		/**
		 * @brief Update all diffusion species.
		 * @param[in] diffusion_reaction_data_i Diffusion reaction data of particle i;
		 * @param[in] dt Time step;
		 */
		void updateDiffusionSpeciesforRungeKutta(DiffusionReactionData& diffusion_reaction_data_i, 
			Real gamma_1, Real gamma_2, Real beta, Real dt)
		{
			for (size_t m = 0; m < species_diffusion_.size(); ++m)
			{
				size_t k = species_diffusion_[m]->diffusion_species_index_;
				Real s_temp = gamma_1 * diffusion_reaction_data_i.species_n_[k] +
					gamma_2 * diffusion_reaction_data_i.species_s_[k] +
					beta * dt * diffusion_reaction_data_i.dspecies_dt_[k];
				diffusion_reaction_data_i.species_n_[k] = s_temp;
			}
		};

		virtual void Update(size_t index_particle_i, Real dt = 0.0) override
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			updateDiffusionSpeciesforRungeKutta(diffusion_reaction_data_i, gamma_1_[RK_step_], gamma_2_[RK_step_], beta_[RK_step_], dt);
		};

	public:
		RelaxationOfAllDifussionSpeciesRungeKutta(BodyType* body)
			: DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>(body) {
			species_diffusion_ = this->material_->getDiffusionSpecies();
		};

		virtual ~RelaxationOfAllDifussionSpeciesRungeKutta() {};

		void setUpRungeKuttaStep(int step) { RK_step_ = step; }
	};

	/**
	* @class RelaxationOfAllReactionsFoward
	* @brief Compute the reaction process of all species by forward splitting
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllReactionsFoward :
		public DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
		/** The reaction model for all reactive species. */
		BaseReactionModel* species_reaction_;
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
			return input * exp(-loss_rate * dt) + production_rate * (1.0 - exp(-loss_rate * dt)) / (loss_rate + 1.0e-30);
		};
		/** Get change rate for all rective species by forward sweeping.
		 * @brief Get change rate for all rective species by backward sweeping.
		 * @param[in] diffusion_reaction_data_i Diffusion Reaction Data.
		 * @param[in] dt Time step size.
		 **/
		void UpdateReactiveSpeciesForward(DiffusionReactionData& diffusion_reaction_data_i, Real dt) {
			IndexVector& reactive_species = species_reaction_->reactive_species_;

			for (size_t m = 0; m != reactive_species.size(); ++m) {
				size_t k = reactive_species[m];
				Real production_rate = species_reaction_->get_production_rates_[k](diffusion_reaction_data_i.species_n_);
				Real loss_rate = species_reaction_->get_loss_rates_[k](diffusion_reaction_data_i.species_n_);
				Real input = diffusion_reaction_data_i.species_n_[k];
				diffusion_reaction_data_i.species_n_[k] = updateAReactionSpecies(input, production_rate, loss_rate, dt);
			}
		};

		virtual void Update(size_t index_particle_i, Real dt = 0.0) override {
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionData& diffusion_reaction_data_i = particles->diffusion_reaction_data_[index_particle_i];
			UpdateReactiveSpeciesForward(diffusion_reaction_data_i, dt);
		};
	public:
		RelaxationOfAllReactionsFoward(BodyType* body)
			: DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>(body) {
			species_reaction_ = this->material_->getReactionModel();
		};
		virtual ~RelaxationOfAllReactionsFoward() {};
	};

	/**
	* @class RelaxationOfAllReactionsBackward
	* @brief Compute the reaction process of all species by backward splitting
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllReactionsBackward :
		public DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
		/** The reaction model for all reactive species. */
		BaseReactionModel* species_reaction_;
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
			return input * exp(-loss_rate * dt) + production_rate * (1.0 - exp(-loss_rate * dt)) / (loss_rate + 1.0e-30);
		};

		/**
		 * @brief Get change rate for all rective species by backward sweeping.
		 * @param[in] diffusion_reaction_data_i Diffusion Reaction Data.
		 * @param[in] dt Time step size.
		 **/
		void UpdateReactiveSpeciesBackward(DiffusionReactionData& diffusion_reaction_data_i, Real dt)
		{
			IndexVector& reactive_species = species_reaction_->reactive_species_;

			for (size_t m = reactive_species.size(); m != 0; --m) {
				size_t k = reactive_species[m - 1];
				Real production_rate = species_reaction_->get_production_rates_[k](diffusion_reaction_data_i.species_n_);
				Real loss_rate = species_reaction_->get_loss_rates_[k](diffusion_reaction_data_i.species_n_);
				Real input = diffusion_reaction_data_i.species_n_[k];
				diffusion_reaction_data_i.species_n_[k] = updateAReactionSpecies(input, production_rate, loss_rate, dt);
			}
		};

		virtual void Update(size_t index_particle_i, Real dt = 0.0) override
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionData& diffusion_reaction_data_i = particles->diffusion_reaction_data_[index_particle_i];
			UpdateReactiveSpeciesBackward(diffusion_reaction_data_i, dt);
		};
	public:
		RelaxationOfAllReactionsBackward(BodyType* body)
			: DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>(body) 
		{
			species_reaction_ = this->material_->getReactionModel();
		};
		virtual ~RelaxationOfAllReactionsBackward() {};
	};
	/**
	 * @class DiffusionBoundaryCondtion
	 * @brief set boudary condition for diffusion problem
	 */
	template <class BodyType, class BaseParticlesType, class BodyPartByParticleType, class BaseMaterialType>
	class DiffusionBoundaryCondtion :
		public DiffusionReactionConstraint<BodyType, BaseParticlesType, BodyPartByParticleType, BaseMaterialType>
	{
	public:
		DiffusionBoundaryCondtion(BodyType* body, BodyPartByParticleType* body_part)
			: DiffusionReactionConstraint<BodyType,
			BaseParticlesType, BodyPartByParticleType, BaseMaterialType>(body, body_part) {}
		virtual ~DiffusionBoundaryCondtion() {};

	};

	/**
	* @class DiffusionBasedMapping
	* @brief Mapping inside of body according to diffusion.
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class DiffusionBasedMapping :
		public DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
	public:
		DiffusionBasedMapping(BodyType* body)
			: DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>(body) {};
	virtual ~DiffusionBasedMapping() {};
	};
}