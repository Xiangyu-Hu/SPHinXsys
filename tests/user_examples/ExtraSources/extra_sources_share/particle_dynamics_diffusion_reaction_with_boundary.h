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
	 * @brief Pure abstract class for initial conditions.
	 */
	template <class DiffusionReactionParticlesType>
	class DiffusionReactionInitialConditionWithBoundary
		: public DiffusionReactionInitialCondition<DiffusionReactionParticlesType>
	{
	public:
		explicit DiffusionReactionInitialConditionWithBoundary(SPHBody& sph_body);
		virtual ~DiffusionReactionInitialConditionWithBoundary() {};

	protected:
		StdLargeVec<Real> heat_flux_;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesNeumannContact
	 * @brief Contact diffusion relaxation with Neumann boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesNeumannContact
		: public RelaxationOfAllDiffusionSpeciesBaseContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_heat_flux_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateNeumannContact(size_t particle_i, size_t particle_j, Real surface_area_ij_Neumann, StdLargeVec<Real>& heat_flux_k);

	public:
		explicit RelaxationOfAllDiffusionSpeciesNeumannContact(ContactRelation& contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesNeumannContact() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H