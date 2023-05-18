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
	class DiffusionReactionInitialConditionWithNeumann
		: public DiffusionReactionInitialCondition<DiffusionReactionParticlesType>
	{
	public:
		explicit DiffusionReactionInitialConditionWithNeumann(SPHBody& sph_body);
		virtual ~DiffusionReactionInitialConditionWithNeumann() {};

	protected:
		StdLargeVec<Real>& heat_flux_;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesNeumannContact
	 * @brief Contact diffusion relaxation with Neumann boundary condition.
	 */
	template <class DiffusionReactionParticlesType, class ContactDiffusionReactionParticlesType>
	class RelaxationOfAllDiffusionSpeciesNeumann
		: public RelaxationOfAllDiffusionSpeciesContact<DiffusionReactionParticlesType, ContactDiffusionReactionParticlesType>
	{
		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_heat_flux_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateNeumannContact(size_t particle_i, size_t particle_j, Real surface_area_ij_Neumann, StdLargeVec<Real>& heat_flux_k);

	public:
		explicit RelaxationOfAllDiffusionSpeciesNeumann(BaseContactRelation& contact_relation);
		virtual ~RelaxationOfAllDiffusionSpeciesNeumann() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H