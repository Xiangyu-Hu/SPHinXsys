/**
 * @file 	particle_dynamics_diffusion_reaction_with_boundary.h
 * @brief 	This is the derived class of particle dynamics with boundary conditions.
 * @author	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_COMPLEXRELATION_H
#define PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_COMPLEXRELATION_H

#include "particle_dynamics_diffusion_reaction.h"

namespace SPH
{
	/**
	 * @class DiffusionReactionInitialCondition
	 * @brief Pure abstract class for initial conditions
	 */
	template <class DiffusionReactionParticlesType>
	class DiffusionReactionInitialConditionWithBoundaryComplexrelation
		: public DiffusionReactionInitialCondition<DiffusionReactionParticlesType>
	{
	public:
		explicit DiffusionReactionInitialConditionWithBoundaryComplexrelation(SPHBody& sph_body);
		virtual ~DiffusionReactionInitialConditionWithBoundaryComplexrelation() {};

	protected:
		StdLargeVec<Real>& heat_flux_;
		StdLargeVec<Real>& convection_;
		StdLargeVec<Real>& T_infinity_;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesSimpleContact
	 * @brief Simple diffusion relaxation process between two contact bodies, which is the base class of three boundary conditions.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesSimpleComplex
		: public LocalDynamics,
		  public DiffusionReactionInnerData<DiffusionReactionParticlesType>,
		  public DiffusionReactionContactData<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
	protected:
		typedef typename DiffusionReactionParticlesType::DiffusionReactionMaterial DiffusionReactionMaterial;
		DiffusionReactionMaterial& material_;
		StdVec<BaseDiffusion*>& all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>*>& gradient_species_;
		StdVec<StdLargeVec<Real>*> diffusion_dt_; //add &

		//void initializeDiffusionChangeRate(size_t particle_i);
		//virtual void updateSpeciesDiffusion(size_t particle_i, Real dt);

	public:
		StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		typedef DiffusionReactionParticlesType InnerParticlesType;
		typedef ComplexRelation BodyRelationType;

		explicit RelaxationOfAllDiffusionSpeciesSimpleComplex(ComplexRelation& complex_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesSimpleComplex() {};
		StdVec<BaseDiffusion*>& AllDiffusions() { return material_.AllDiffusions(); };

		//virtual void interaction(size_t index_i, Real dt = 0.0) = 0;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesDirichletContact
	 * @brief Contact diffusion relaxation with Dirichlet boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesDirichletComplex
		: public RelaxationOfAllDiffusionSpeciesSimpleComplex<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

	protected:
		void getDiffusionChangeRateDirichletComplex(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij,
			const StdVec<StdLargeVec<Real>*>& gradient_species_k);

	public:
		//typedef ContactRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesDirichletComplex(ComplexRelation& complex_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesDirichletComplex() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesNeumannContact
	 * @brief Contact diffusion relaxation with Neumann boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesNeumannComplex
		: public RelaxationOfAllDiffusionSpeciesSimpleComplex<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec<StdLargeVec<Real>*> contact_heat_flux_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateNeumannComplex(size_t particle_i, size_t particle_j, Real surface_area_ij_Neumann, StdLargeVec<Real>& heat_flux_k);

	public:
		//typedef ContactRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesNeumannComplex(ComplexRelation& complex_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesNeumannComplex() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRobinContact
	 * @brief Contact diffusion relaxation with Robin boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesRobinComplex
		: public RelaxationOfAllDiffusionSpeciesSimpleComplex<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;

		StdVec<StdLargeVec<Real>*> contact_convection_;
		StdVec<StdLargeVec<Real>*> contact_T_infinity_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateRobinComplex(size_t particle_i, size_t particle_j, Real surface_area_ij_Robin, StdLargeVec<Real>& convection_k, StdLargeVec<Real>& T_infinity_k);

	public:
		//typedef ContactRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesRobinComplex(ComplexRelation& complex_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesRobinComplex() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class ComplexInteraction
	 * @brief
	 */
	template <typename... DiffusionRelaxationType>
	class ComplexInteraction;

	/*template <>
	class ComplexInteraction<> : public LocalDynamics
	{
	public:
		template <class BodyRelationType>
		ComplexInteraction(BodyRelationType& body_relation)
			: LocalDynamics(body_relation.getDynamicsRange()) {};

		void interaction(size_t index_i, Real dt = 0.0) {};
	};*/

	template <>
	class ComplexInteraction<>
	{
	public:
		template <class BodyRelation>
		ComplexInteraction(BodyRelation& body_relation) {};

		void interaction(size_t index_i, Real dt = 0.0) {};
	};

	template <class DiffusionRelaxationFirst, class... DiffusionRelaxationOthers>
	class ComplexInteraction<DiffusionRelaxationFirst, DiffusionRelaxationOthers...> : public DiffusionRelaxationFirst
	{
	protected:
		//DiffusionRelaxationFirst first_diffusion_relaxation_;
		ComplexInteraction<DiffusionRelaxationOthers...> others_diffusion_relaxation_;

	public:
		template <typename... ComplexArgs>
		ComplexInteraction(BaseInnerRelation& inner_relation, ComplexArgs &&...complex_args)
			: DiffusionRelaxationFirst(inner_relation),
			others_diffusion_relaxation_(std::forward<ComplexArgs>(complex_args)...) {};

		typedef typename DiffusionRelaxationFirst::InnerParticlesType InnerParticlesTypeComplex;
		typedef typename DiffusionRelaxationFirst::BodyRelationType BodyRelationTypeComplex;
		typedef typename DiffusionRelaxationOthers::BodyRelationType BodyRelationTypeComplex;

		void interaction(size_t index_i, Real dt = 0.0)
		{
			DiffusionRelaxationFirst::interaction(index_i, dt);
			others_diffusion_relaxation_.interaction(index_i, dt);
		};
	};

	/**
	 * @class InitializationRKComplex
	 * @brief initialization of a runge-kutta integration scheme
	 */
	template <class DiffusionReactionParticlesType>
	class InitializationRKComplex : public LocalDynamics,
		public DiffusionReactionSimpleData<DiffusionReactionParticlesType>
	{
	protected:
		typename DiffusionReactionParticlesType::DiffusionReactionMaterial& material_;
		StdVec<BaseDiffusion*>& all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>>& diffusion_species_s_;

	public:
		InitializationRKComplex(StdVec<StdLargeVec<Real>>& diffusion_species_s, SPHBody& sph_body);
		virtual ~InitializationRKComplex() {};

		void update(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class SecondStageRK2Complex
	 * @brief the second stage of the 2nd-order Runge-Kutta scheme
	 */
	template <class FirstStageType>
	class SecondStageRK2Complex : public FirstStageType
	{
	protected:
		StdVec<StdLargeVec<Real>>& diffusion_species_s_;
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;

	public:
		template <typename... ComplexArgsType>
		SecondStageRK2Complex(StdVec<StdLargeVec<Real>>& diffusion_species_s, typename FirstStageType::BodyRelationTypeComplex& body_relation, ComplexArgsType &&... complex_agrs)
			: FirstStageType(body_relation, std::forward<ComplexArgsType>(complex_agrs)...), diffusion_species_s_(diffusion_species_s) {};
		virtual ~SecondStageRK2Complex() {};
	};

	template <class FirstStageType>
	class RelaxationOfAllDiffusionSpeciesRK2Complex : public BaseDynamics<void>
	{
	protected:
		/** Intermediate Value */
		StdVec<StdLargeVec<Real>> diffusion_species_s_;
		SimpleDynamics<InitializationRKComplex<typename FirstStageType::InnerParticlesTypeComplex>> rk2_initialization_;
		InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
		InteractionWithUpdate<SecondStageRK2Complex<FirstStageType>> rk2_2nd_stage_;
		StdVec<BaseDiffusion*> all_diffusions_;

	public:
		template <typename... ComplexArgsType>
		explicit RelaxationOfAllDiffusionSpeciesRK2Complex(typename FirstStageType::BodyRelationTypeComplex& body_relation, ComplexArgsType &&... complex_agrs)
			: BaseDynamics<void>(body_relation.getSPHBody()),
			rk2_initialization_(diffusion_species_s_, body_relation.getSPHBody()),
			rk2_1st_stage_(body_relation, std::forward<ComplexArgsType>(complex_agrs)...),
			rk2_2nd_stage_(diffusion_species_s_, body_relation, std::forward<ComplexArgsType>(complex_agrs)...),
			all_diffusions_(rk2_1st_stage_.AllDiffusions())
		{
			diffusion_species_s_.resize(all_diffusions_.size());
			StdVec<std::string>& all_species_names = rk2_1st_stage_.getParticles()->AllSpeciesNames();
			for (size_t i = 0; i != all_diffusions_.size(); ++i)
			{
				// register diffusion species intermediate
				size_t diffusion_species_index = all_diffusions_[i]->diffusion_species_index_;
				std::string& diffusion_species_name = all_species_names[diffusion_species_index];
				rk2_1st_stage_.getParticles()->registerVariable(diffusion_species_s_[i], diffusion_species_name + "Intermediate");
			}
		}
		virtual ~RelaxationOfAllDiffusionSpeciesRK2Complex() {};

		virtual void exec(Real dt = 0.0) override;
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
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H