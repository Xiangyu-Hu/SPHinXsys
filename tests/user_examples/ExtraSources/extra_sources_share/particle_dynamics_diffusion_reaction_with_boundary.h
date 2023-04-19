/**
 * @file 	particle_dynamics_diffusion_reaction_with_boundary.h
 * @brief 	This is the derived class of particle dynamics with boundary conditions.
 * @author	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H
#define PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H

#include "particle_dynamics_diffusion_reaction.h"

namespace SPH
{

	/**
	 * @class DiffusionReactionInitialCondition
	 * @brief Pure abstract class for initial conditions
	 */
	template <class DiffusionReactionParticlesType>
	class DiffusionReactionInitialConditionWithBoundary
		: public DiffusionReactionInitialCondition<DiffusionReactionParticlesType>
	{
	public:
		explicit DiffusionReactionInitialConditionWithBoundary(SPHBody& sph_body);
		virtual ~DiffusionReactionInitialConditionWithBoundary() {};

	protected:
		StdLargeVec<Real>& heat_flux_;
		StdLargeVec<Real>& convection_;
		StdLargeVec<Real>& T_infinity_;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesWithDirichlet
	 * @brief Contact diffusion relaxation with Dirichlet boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesWithDirichlet
		: public RelaxationOfAllDiffusionSpeciesInner<DiffusionReactionParticlesType>,
		public DiffusionReactionContactData<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
	protected:
		void getDiffusionChangeRateWithDirichlet(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij,
			const StdVec<StdLargeVec<Real>*>& gradient_species_k);

	public:
		StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		typedef ComplexRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesWithDirichlet(ComplexRelation& complex_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesWithDirichlet() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesWithNeumann
	 * @brief Contact diffusion relaxation with Neumann boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesWithNeumann
		: public RelaxationOfAllDiffusionSpeciesWithDirichlet<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec<StdLargeVec<Real>*> contact_heat_flux_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateWithNeumann(size_t particle_i, size_t particle_j, Real surface_area_ij_Neumann, StdLargeVec<Real>& heat_flux_k);

	public:
		//StdVec<StdLargeVec<Real>*>& GradientSpecies() { return this->contact_gradient_species_; };

		typedef ComplexRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesWithNeumann(ComplexRelation& complex_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesWithNeumann() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesWithRobin
	 * @brief Contact diffusion relaxation with Robin boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesWithRobin
		: public RelaxationOfAllDiffusionSpeciesWithDirichlet<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		
		StdVec<StdLargeVec<Real>*> contact_convection_;
		StdVec<StdLargeVec<Real>*> contact_T_infinity_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateWithRobin(size_t particle_i, size_t particle_j, Real surface_area_ij_Robin, StdLargeVec<Real>& convection_k, StdLargeVec<Real>& T_infinity_k);

	public:
		//StdVec<StdLargeVec<Real>*>& GradientSpecies() { return this->contact_gradient_species_; };

		typedef ComplexRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesWithRobin(ComplexRelation& complex_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesWithRobin() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class UpdateUnitVectorNormalToBoundary
	 * @brief Update the normal vector of surface.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class UpdateUnitVectorNormalToBoundary
		: public LocalDynamics,
		public DiffusionReactionInnerData<DiffusionReactionParticlesType>,
		public DiffusionReactionContactData<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
	protected:
		StdLargeVec<Vecd>& normal_vector_;

	public:
		UpdateUnitVectorNormalToBoundary(ComplexRelation& body_complex_relation);
		virtual ~UpdateUnitVectorNormalToBoundary() {};
		void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class InitializationRKWithBoundary
	 * @brief initialization of a runge-kutta integration scheme
	 */
	template <class DiffusionReactionParticlesType>
	class InitializationRKWithBoundary : public LocalDynamics,
		public DiffusionReactionSimpleData<DiffusionReactionParticlesType>
	{
	protected:
		typename DiffusionReactionParticlesType::DiffusionReactionMaterial& material_;
		StdVec<BaseDiffusion*>& all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_s_;

	public:
		InitializationRKWithBoundary(SPHBody& sph_body, StdVec<StdLargeVec<Real>*>& diffusion_species_s);
		virtual ~InitializationRKWithBoundary() {};

		void update(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class SecondStageRK2WithBoundary
	 * @brief the second stage of the 2nd-order Runge-Kutta scheme
	 */
	template <class FirstStageType>
	class SecondStageRK2WithBoundary : public FirstStageType
	{
	protected:
		StdVec<StdLargeVec<Real>*>& diffusion_species_s_;
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;

	public:
		SecondStageRK2WithBoundary(typename FirstStageType::BodyRelationType& body_relation,
			StdVec<StdLargeVec<Real>*>& diffusion_species_s);
		virtual ~SecondStageRK2WithBoundary() {};
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRK2WithBoundary
	 * @brief Compute the diffusion relaxation process of all species
	 * with second order Runge-Kutta time stepping
	 */
	template <class FirstStageType>
	class RelaxationOfAllDiffusionSpeciesRK2WithBoundary : public BaseDynamics<void>
	{
	protected:
		/** Intermediate Value */

		StdVec<StdLargeVec<Real>*> diffusion_species_s_;
		SimpleDynamics<InitializationRKWithBoundary<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
		InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
		InteractionWithUpdate<SecondStageRK2WithBoundary<FirstStageType>> rk2_2nd_stage_;
		StdVec<BaseDiffusion*> all_diffusions_;

	public:
		explicit RelaxationOfAllDiffusionSpeciesRK2WithBoundary(typename FirstStageType::BodyRelationType& body_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesRK2WithBoundary() {};

		virtual void exec(Real dt = 0.0) override;
	};



	/**
	 * @class InitializationRK02
	 * @brief initialization of a runge-kutta integration scheme
	 */
	template <class DiffusionReactionParticlesType>
	class InitializationRK02 : public LocalDynamics,
		public DiffusionReactionSimpleData<DiffusionReactionParticlesType>
	{
	protected:
		typename DiffusionReactionParticlesType::DiffusionReactionMaterial& material_;
		StdVec<BaseDiffusion*>& all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>>& diffusion_species_s_;

	public:
		InitializationRK02(SPHBody& sph_body, StdVec<StdLargeVec<Real>>& diffusion_species_s);
		virtual ~InitializationRK02() {};

		void update(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class SecondStageRK202
	 * @brief the second stage of the 2nd-order Runge-Kutta scheme
	 */
	template <class FirstStageType>
	class SecondStageRK202 : public FirstStageType
	{
	protected:
		StdVec<StdLargeVec<Real>>& diffusion_species_s_;
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;

	public:
		SecondStageRK202(typename FirstStageType::BodyRelationType& body_relation,
			StdVec<StdLargeVec<Real>>& diffusion_species_s);
		virtual ~SecondStageRK202() {};
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRK2
	 * @brief Compute the diffusion relaxation process of all species
	 * with second order Runge-Kutta time stepping
	 */
	template <class FirstStageType>
	class RelaxationOfAllDiffusionSpeciesRK202 : public BaseDynamics<void>
	{
	protected:
		/** Intermediate Value */
		StdVec<StdLargeVec<Real>> diffusion_species_s_;
		SimpleDynamics<InitializationRK02<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
		InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
		InteractionWithUpdate<SecondStageRK202<FirstStageType>> rk2_2nd_stage_;
		StdVec<BaseDiffusion*> all_diffusions_;

	public:
		explicit RelaxationOfAllDiffusionSpeciesRK202(typename FirstStageType::BodyRelationType& body_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesRK202() {};

		virtual void exec(Real dt = 0.0) override;
	};
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H