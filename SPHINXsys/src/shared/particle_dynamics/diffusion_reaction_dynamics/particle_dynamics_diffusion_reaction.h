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
*/

#pragma once

#include "all_particle_dynamics.h"
#include "diffusion_reaction_particles.h"
#include "diffusion_reaction.h"


namespace SPH
{
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionBase = ParticleDynamics<Real, BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionReactionSimple = ParticleDynamicsSimple<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionInner = ParticleDynamicsInner<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseMaterialType>>;

	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	using DiffusionInnerWithUpdate = ParticleDynamicsInnerWithUpdate<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, DiffusionReactionMaterial<BaseMaterialType>>;

	template <class BodyType, class BaseParticlesType, class BodyPartByParticleType, class BaseMaterialType>
	using DiffusionBoundaryConditionConstraint = ConstraintByParticle<BodyType,
		DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>, BodyPartByParticleType, DiffusionReactionMaterial<BaseMaterialType>>;

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
	protected:
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;
			Neighborhood& neighborhood = (*this->reference_inner_configuration_)[index_particle_i];

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			StdLargeVec<BaseParticleData>& base_particle_data = particles->base_particle_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];

			material->initializeDiffusionChangeRate(diffusion_reaction_data_i);
			NeighborList& neighors = std::get<0>(neighborhood);
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle* neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real Vol_j = base_particle_data[index_particle_j].Vol_;
				DiffusionReactionData& diffusion_reaction_data_j = diffusion_reaction_data[index_particle_j];
	
				const Vecd& gradi_ij = particles->getKernelGradient(index_particle_i, index_particle_j, neighboring_particle->dW_ij_, e_ij);
				Real area_ij = 2.0 * base_particle_data[index_particle_j].Vol_ * dot(gradi_ij, e_ij) / neighboring_particle->r_ij_;
				material->getDiffusionChangeRate(index_particle_i, index_particle_j, e_ij, area_ij,
					diffusion_reaction_data_i, diffusion_reaction_data_j);
			}
		};
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			material->updateDiffusionSpecies(diffusion_reaction_data_i, dt);
		};
	public:
		RelaxationOfAllDifussionSpecies(BodyType* body)
			: DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>(body) {};
		virtual ~RelaxationOfAllDifussionSpecies() {};
	};

	/**
	* @class RungeKuttaRelaxationOfAllDifussionSpecies
	* @brief Compute the diffusion relaxation process of all species using 4th Runge-Kuttas cheme
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class TotalLagrangianRungeKuttaRelaxationOfAllDifussionSpecies
		: public DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>
	{
	private: 
		Real gamma_1_[5] = {0.0, 0.0, 0.12109848,-3.8438337, 0.5463709};
		Real gamma_2_[5] = {0.0, 1.0, 0.7217817,  2.1212093, 0.1986530};
		Real beta_[5] 	 = {0.0, 1.1937439, 0.0992799, 1.1316780, 0.3106657};
		Real delta_[4] 	 = {1.0, 0.2176833, 1.0658413, 0.0};
		int RK_step_;
	protected:
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;
			Neighborhood& neighborhood = (*this->reference_inner_configuration_)[index_particle_i];

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			StdLargeVec<BaseParticleData>& base_particle_data = particles->base_particle_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			if(RK_step_ == 1)
			{
				material->initializeStageforRungeKutta(diffusion_reaction_data_i);
			}
			material->updateStageforRungeKutta(diffusion_reaction_data_i, delta_[RK_step_-1]);
			material->initializeDiffusionChangeRate(diffusion_reaction_data_i);
			NeighborList& neighors = std::get<0>(neighborhood);
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle* neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real Vol_j = base_particle_data[index_particle_j].Vol_;
				DiffusionReactionData& diffusion_reaction_data_j = diffusion_reaction_data[index_particle_j];
	
				const Vecd& gradi_ij = particles->getKernelGradient(index_particle_i, index_particle_j, neighboring_particle->dW_ij_, e_ij);
				Real area_ij = 2.0 * base_particle_data[index_particle_j].Vol_ * dot(gradi_ij, e_ij) / neighboring_particle->r_ij_;
				material->getDiffusionChangeRate(index_particle_i, index_particle_j, e_ij, area_ij,
					diffusion_reaction_data_i, diffusion_reaction_data_j);
			}
		};

		virtual void Update(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			material->updateDiffusionSpeciesforRungeKutta(diffusion_reaction_data_i, gamma_1_[RK_step_], gamma_2_[RK_step_], beta_[RK_step_], dt);
		};

	public:
		TotalLagrangianRungeKuttaRelaxationOfAllDifussionSpecies(BodyType* body)
			: DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>(body){};
		virtual ~TotalLagrangianRungeKuttaRelaxationOfAllDifussionSpecies() {};

		void setUpRungeKuttaStep(int step){RK_step_ = step;}
	};

	/**
	* @class RungeKuttaRelaxationOfAllDifussionSpecies
	* @brief Compute the diffusion relaxation process of all species using 4th Runge-Kuttas cheme
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class UpdateLagrangianRungeKuttaRelaxationOfAllDifussionSpecies
		: public DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>
	{
	private: 
		Real gamma_1_[5] = {0.0, 0.0, 0.12109848,-3.8438337, 0.5463709};
		Real gamma_2_[5] = {0.0, 1.0, 0.7217817,  2.1212093, 0.1986530};
		Real beta_[5] 	 = {0.0, 1.1937439, 0.0992799, 1.1316780, 0.3106657};
		Real delta_[4] 	 = {1.0, 0.2176833, 1.0658413, 0.0};
		int RK_step_;
	protected:
		virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;
			Neighborhood& neighborhood = (*this->current_inner_configuration_)[index_particle_i];

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			StdLargeVec<BaseParticleData>& base_particle_data = particles->base_particle_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			if(RK_step_ == 1)
			{
				material->initializeStageforRungeKutta(diffusion_reaction_data_i);
			}
			material->updateStageforRungeKutta(diffusion_reaction_data_i, delta_[RK_step_-1]);
			material->initializeDiffusionChangeRate(diffusion_reaction_data_i);
			NeighborList& neighors = std::get<0>(neighborhood);
			for (size_t n = 0; n < neighors.size(); ++n)
			{
				NeighboringParticle* neighboring_particle = neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				Vecd& e_ij = neighboring_particle->e_ij_;
				Real Vol_j = base_particle_data[index_particle_j].Vol_;
				DiffusionReactionData& diffusion_reaction_data_j = diffusion_reaction_data[index_particle_j];
	
				const Vecd& gradi_ij = particles->getKernelGradient(index_particle_i, index_particle_j, neighboring_particle->dW_ij_, e_ij);
				Real area_ij = 2.0 * base_particle_data[index_particle_j].Vol_ * dot(gradi_ij, e_ij) / neighboring_particle->r_ij_;
				material->getDiffusionChangeRate(index_particle_i, index_particle_j, e_ij, area_ij,
					diffusion_reaction_data_i, diffusion_reaction_data_j);
			}
		};

		virtual void Update(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;

			StdLargeVec<DiffusionReactionData>& diffusion_reaction_data = particles->diffusion_reaction_data_;
			DiffusionReactionData& diffusion_reaction_data_i = diffusion_reaction_data[index_particle_i];
			material->updateDiffusionSpeciesforRungeKutta(diffusion_reaction_data_i, gamma_1_[RK_step_], gamma_2_[RK_step_], beta_[RK_step_], dt);
		};

	public:
		UpdateLagrangianRungeKuttaRelaxationOfAllDifussionSpecies(BodyType* body)
			: DiffusionInnerWithUpdate<BodyType, BaseParticlesType, BaseMaterialType>(body){};
		virtual ~UpdateLagrangianRungeKuttaRelaxationOfAllDifussionSpecies() {};

		void setUpRungeKuttaStep(int step){RK_step_ = step;}
	};
	/**
	* @class RelaxationOfAllReactionsFoward
	* @brief Compute the reaction process of all species by forward splitting
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RelaxationOfAllReactionsFoward :
		public DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override {
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;
	
			DiffusionReactionData& diffusion_reaction_data_i = particles->diffusion_reaction_data_[index_particle_i];
			material->UpdateReactiveSpeciesForward(diffusion_reaction_data_i, dt);
		};
	public:
		RelaxationOfAllReactionsFoward(BodyType* body)
			: DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>(body) {};
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
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override 
		{
			DiffusionReactionParticles<BaseParticlesType, BaseMaterialType>* particles = this->particles_;
			DiffusionReactionMaterial<BaseMaterialType>* material = this->material_;

			DiffusionReactionData& diffusion_reaction_data_i = particles->diffusion_reaction_data_[index_particle_i];
			material->UpdateReactiveSpeciesBackward(diffusion_reaction_data_i, dt);
		};
	public:
		RelaxationOfAllReactionsBackward(BodyType* body)
			: DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>(body) {};
		virtual ~RelaxationOfAllReactionsBackward() {};
	};
	/**
	 * @class SetDiffusionBoundaryCondtion
	 * @brief set boudary condition for diffusion problem
	 */
	template <class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType>
	class SetDiffusionBoundaryCondtion :
		public DiffusionBoundaryConditionConstraint<BodyType,ParticlesType, BodyPartByParticleType, MaterialType>
	{
	protected:
		virtual void ConstraintAParticle(size_t index_particle_i, Real dt = 0.0) override {};
	public:
		SetDiffusionBoundaryCondtion(BodyType* body, BodyPartByParticleType* body_part)
			: DiffusionBoundaryConditionConstraint<BodyType,ParticlesType, BodyPartByParticleType, MaterialType>(body, body_part) {}
		virtual ~SetDiffusionBoundaryCondtion() {};

	};

	/**
	 * @class ConstraintDiffusionBoundaryCondtion
	 * @brief constrain boudary condition for diffusion problem
	 */
	template <class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType>
	class ConstraintDiffusionBoundaryCondtion :
		public DiffusionBoundaryConditionConstraint<BodyType,ParticlesType, BodyPartByParticleType, MaterialType>
	{
	protected:
		virtual void ConstraintAParticle(size_t index_particle_i, Real dt = 0.0) override {};
	public:
		ConstraintDiffusionBoundaryCondtion(BodyType* body, BodyPartByParticleType* body_part)
			: DiffusionBoundaryConditionConstraint<BodyType,ParticlesType, BodyPartByParticleType, MaterialType>(body, body_part) {}
		virtual ~ConstraintDiffusionBoundaryCondtion() {};

	};

	/**
	* @class computeFiberandSheetDirection
	* @brief Compute the fiber and sheet direction by using the coupling of diffusion tensor and
	*		 Level-set method
	*		The diffusion tensor will determine the rotation angle and Levelset will determine
	*		the face norm.
	*/
	template<class BodyType, class BaseParticlesType, class BaseMaterialType>
	class RuleBasedFibreandSheetConstruction  :
		public DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>
	{
	protected:
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override{ };
	public:
		RuleBasedFibreandSheetConstruction(BodyType* body)
			: DiffusionReactionSimple<BodyType, BaseParticlesType, BaseMaterialType>(body) { };
	virtual ~RuleBasedFibreandSheetConstruction() {};
	};
}
