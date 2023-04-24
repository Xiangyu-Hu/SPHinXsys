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
	 * @class RelaxationOfAllDiffusionSpeciesSimpleContact
	 * @brief Simple diffusion relaxation process between two contact bodies, which is the base class of three boundary conditions.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesSimpleContact
		: public LocalDynamics,
		public DiffusionReactionContactData<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
	protected:
		typedef typename DiffusionReactionParticlesType::DiffusionReactionMaterial DiffusionReactionMaterial;
		DiffusionReactionMaterial& material_;
		StdVec<BaseDiffusion*>& all_diffusions_;
		StdVec<StdLargeVec<Real>*>& diffusion_species_;
		StdVec<StdLargeVec<Real>*>& gradient_species_;
		StdVec<StdLargeVec<Real>*> diffusion_dt_;

		void initializeDiffusionChangeRate(size_t particle_i);
		//virtual void updateSpeciesDiffusion(size_t particle_i, Real dt);

	public:
		StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		typedef DiffusionReactionParticlesType InnerParticlesType;
		typedef ContactRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesSimpleContact(ContactRelation& contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesSimpleContact() {};
		StdVec<BaseDiffusion*>& AllDiffusions() { return material_.AllDiffusions(); };

		virtual void interaction(size_t index_i, Real dt = 0.0) = 0;  // right or not?
	};


	/**
	 * @class RelaxationOfAllDiffusionSpeciesDirichletContact
	 * @brief Contact diffusion relaxation with Dirichlet boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesDirichletContact
		: public RelaxationOfAllDiffusionSpeciesSimpleContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

	protected:
		void getDiffusionChangeRateDirichletContact(size_t particle_i, size_t particle_j, Vecd& e_ij, Real surface_area_ij,
			const StdVec<StdLargeVec<Real>*>& gradient_species_k);

	public:
		//typedef ContactRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesDirichletContact(ContactRelation& contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesDirichletContact() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesNeumannContact
	 * @brief Contact diffusion relaxation with Neumann boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesNeumannContact
		: public RelaxationOfAllDiffusionSpeciesSimpleContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec<StdLargeVec<Real>*> contact_heat_flux_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateNeumannContact(size_t particle_i, size_t particle_j, Real surface_area_ij_Neumann, StdLargeVec<Real>& heat_flux_k);

	public:
		//typedef ContactRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesNeumannContact(ContactRelation& contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesNeumannContact() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRobinContact
	 * @brief Contact diffusion relaxation with Robin boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesRobinContact
		: public RelaxationOfAllDiffusionSpeciesSimpleContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		//StdVec<StdVec<StdLargeVec<Real>*>> contact_gradient_species_;

		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_Vol_;

		StdVec<StdLargeVec<Real>*> contact_convection_;
		StdVec<StdLargeVec<Real>*> contact_T_infinity_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateRobinContact(size_t particle_i, size_t particle_j, Real surface_area_ij_Robin, StdLargeVec<Real>& convection_k, StdLargeVec<Real>& T_infinity_k);

	public:
		//typedef ContactRelation BodyRelationType;
		explicit RelaxationOfAllDiffusionSpeciesRobinContact(ContactRelation& contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesRobinContact() {};

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
		template <class BodyRelationType>
		ComplexInteraction(BodyRelationType& body_relation) {};

		void interaction(size_t index_i, Real dt = 0.0) {};
	};


	template <class DiffusionRelaxationFirst, class... DiffusionRelaxationOthers>
	class ComplexInteraction<DiffusionRelaxationFirst, DiffusionRelaxationOthers...> : public DiffusionRelaxationFirst
	{
	protected:
		//DiffusionRelaxationFirst first_diffusion_relaxation_;
		ComplexInteraction<DiffusionRelaxationOthers...> others_diffusion_relaxation_;

	public:
		template <typename... ContactArgs>
		ComplexInteraction(BaseInnerRelation& inner_relation, BaseContactRelation& contact_relation, ContactArgs &&...contact_args)
			: DiffusionRelaxationFirst(inner_relation),
			others_diffusion_relaxation_(contact_relation, std::forward<ContactArgs>(contact_args)...) {};

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
		SecondStageRK2Complex(StdVec<StdLargeVec<Real>>& diffusion_species_s, typename FirstStageType::BodyRelationType& body_relation);
		virtual ~SecondStageRK2Complex() {};

		void update(size_t index_i, Real dt = 0.0); // add by Zhao
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRK2Complex
	 * @brief Compute the diffusion relaxation process of all species
	 * with second order Runge-Kutta time stepping
	 */
	//template <class FirstStageType>
	//class RelaxationOfAllDiffusionSpeciesRK2Complex : public BaseDynamics<void>
	//{
	//protected:
	//	/** Intermediate Value */
	//	StdVec<StdLargeVec<Real>> diffusion_species_s_;
	//	SimpleDynamics<InitializationRKComplex<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
	//	InteractionWithUpdate<FirstStageType> rk2_1st_stage_;
	//	InteractionWithUpdate<SecondStageRK2Complex<FirstStageType>> rk2_2nd_stage_;
	//	StdVec<BaseDiffusion*> all_diffusions_;

	//public:
	//	template <typename... ContactArgsType>
	//	explicit RelaxationOfAllDiffusionSpeciesRK2Complex(typename FirstStageType::BodyRelationType& body_relation, ContactArgsType &&... contact_agrs)
	//		: BaseDynamics<void>(body_relation.getSPHBody()),
	//		rk2_initialization_(body_relation.getSPHBody(), diffusion_species_s_),
	//		rk2_1st_stage_(body_relation),
	//		rk2_2nd_stage_(body_relation, diffusion_species_s_),
	//		all_diffusions_(rk2_1st_stage_.AllDiffusions())
	//	{
	//		diffusion_species_s_.resize(all_diffusions_.size());
	//		StdVec<std::string>& all_species_names = rk2_1st_stage_.getParticles()->AllSpeciesNames();
	//		for (size_t i = 0; i != all_diffusions_.size(); ++i)
	//		{
	//			// register diffusion species intermediate
	//			size_t diffusion_species_index = all_diffusions_[i]->diffusion_species_index_;
	//			std::string& diffusion_species_name = all_species_names[diffusion_species_index];
	//			rk2_1st_stage_.getParticles()->registerVariable(diffusion_species_s_[i], diffusion_species_name + "Intermediate");
	//		}
	//	}
	//	virtual ~RelaxationOfAllDiffusionSpeciesRK2Complex() {};

	//	virtual void exec(Real dt = 0.0) override;
	//};

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
			rk2_initialization_(diffusion_species_s_, body_relation.getSPHBody(), std::forward<ContactArgsType>(contact_agrs)...),
			rk2_1st_stage_(body_relation, std::forward<ContactArgsType>(contact_agrs)...),
			rk2_2nd_stage_(diffusion_species_s_, body_relation, std::forward<ContactArgsType>(contact_agrs)...),
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