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
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H