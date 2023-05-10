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
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	using DiffusionReactionContactDataWithBoundary =
		DataDelegateContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>;

	/**
	 * @class RelaxationOfAllDiffusionSpeciesSimpleContact
	 * @brief Simple diffusion relaxation process between two contact bodies, which is the base class of three boundary conditions.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesBaseContact
		: public LocalDynamics,
		  public DiffusionReactionContactDataWithBoundary<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
	protected:
		typedef typename DiffusionReactionParticlesType::DiffusionReactionMaterial DiffusionReactionMaterial;
		DiffusionReactionMaterial& material_;
		StdVec<BaseDiffusion*>& all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>*>& gradient_species_;
		StdVec<StdLargeVec<Real>*> diffusion_dt_;

	public:
		StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		typedef DiffusionReactionParticlesType InnerParticlesType;
		typedef BaseContactRelation BodyRelationType;

		explicit RelaxationOfAllDiffusionSpeciesBaseContact(BaseContactRelation &contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesBaseContact() {};
		StdVec<BaseDiffusion*>& AllDiffusions() { return material_.AllDiffusions(); };

		virtual void interaction(size_t index_i, Real dt = 0.0) = 0;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesDirichletContact
	 * @brief Contact diffusion relaxation with Dirichlet boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesDirichletContact
		: public RelaxationOfAllDiffusionSpeciesBaseContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
	protected:
		void getDiffusionChangeRateDirichletContact(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij,
			const StdVec<StdLargeVec<Real>*>& gradient_species_k);

	public:
		explicit RelaxationOfAllDiffusionSpeciesDirichletContact(ContactRelation& contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesDirichletContact() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class ComplexInteraction
	 * @brief A class that integrates multiple boundary conditions.
	 */
	template <typename... InteractionType>
	class ComplexInteraction;

	template <>
	class ComplexInteraction<>
	{
	public:
		ComplexInteraction() {};

		void interaction(size_t index_i, Real dt = 0.0) {};
	};

	template <class FirstInteraction, class... OtherInteractions>
    class ComplexInteraction<FirstInteraction, OtherInteractions...> : public FirstInteraction
	{
	protected:
		ComplexInteraction<OtherInteractions...> other_interaction_;

	public:

		template <class FirstRelationType, typename... OtherRelationTypes>
		explicit ComplexInteraction(FirstRelationType& body_relation, OtherRelationTypes &&...other_relations)
			: FirstInteraction(body_relation),
			other_interaction_(std::forward<OtherRelationTypes>(other_relations)...) {};

		void interaction(size_t index_i, Real dt = 0.0)
		{
			FirstInteraction::interaction(index_i, dt);
			other_interaction_.interaction(index_i, dt);
		};
	};

	/**
	 * @class InitializationRKComplex
	 * @brief Initialization of a runge-kutta integration scheme.
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
		InitializationRKComplex(SPHBody& sph_body, StdVec<StdLargeVec<Real>>& diffusion_species_s);
		virtual ~InitializationRKComplex() {};

		void update(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class SecondStageRK2Complex
	 * @brief The second stage of the 2nd-order Runge-Kutta scheme.
	 */
	template <class FirstStageType>
	class SecondStageRK2Complex : public FirstStageType
	{
	protected:
		StdVec<StdLargeVec<Real>>& diffusion_species_s_;
		virtual void updateSpeciesDiffusion(size_t particle_i, Real dt) override;

	public:
		template <typename... ContactArgsType>
		SecondStageRK2Complex(typename FirstStageType::BodyRelationType& body_relation, StdVec<StdLargeVec<Real>>& diffusion_species_s, ContactArgsType &&... contact_agrs)
			: FirstStageType(body_relation, std::forward<ContactArgsType>(contact_agrs)...),
			diffusion_species_s_(diffusion_species_s) {};
		virtual ~SecondStageRK2Complex() {};
	};

	template <class FirstStageType>
	class RelaxationOfAllDiffusionSpeciesRK2Complex : public BaseDynamics<void>
	{
	protected:
		/** Intermediate Value */
		StdVec<StdLargeVec<Real>> diffusion_species_s_;
		SimpleDynamics<InitializationRKComplex<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
		InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
		InteractionWithUpdate<SecondStageRK2Complex<FirstStageType>> rk2_2nd_stage_;
		StdVec<BaseDiffusion*> all_diffusions_;

	public:
		template <typename... ContactArgsType>
		explicit RelaxationOfAllDiffusionSpeciesRK2Complex(typename FirstStageType::BodyRelationType& body_relation, ContactArgsType &&... contact_agrs)
			: BaseDynamics<void>(body_relation.getSPHBody()),
			rk2_initialization_(body_relation.getSPHBody(), diffusion_species_s_),
			rk2_1st_stage_(body_relation, std::forward<ContactArgsType>(contact_agrs)...),
			rk2_2nd_stage_(body_relation, diffusion_species_s_, std::forward<ContactArgsType>(contact_agrs)...),
			all_diffusions_(rk2_1st_stage_.AllDiffusions())
		{
			diffusion_species_s_.resize(all_diffusions_.size());
			StdVec<std::string>& all_species_names = rk2_1st_stage_.getParticles()->AllSpeciesNames();
			for (size_t i = 0; i != all_diffusions_.size(); ++i)
			{
				// Register diffusion species intermediate
				size_t diffusion_species_index = all_diffusions_[i]->diffusion_species_index_;
				std::string& diffusion_species_name = all_species_names[diffusion_species_index];
				rk2_1st_stage_.getParticles()->registerVariable(diffusion_species_s_[i], diffusion_species_name + "Intermediate");
			}
		}
		virtual ~RelaxationOfAllDiffusionSpeciesRK2Complex() {};

		virtual void exec(Real dt = 0.0) override;
	};
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H