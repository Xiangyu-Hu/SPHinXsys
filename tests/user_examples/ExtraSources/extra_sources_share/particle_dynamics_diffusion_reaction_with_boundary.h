/**
 * @file 	particle_dynamics_diffusion_reaction_with_boundary.h
 * @brief 	This is the derived class of particle dynamics with boundary conditions.
 * @author	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H
#define PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H

#include "diffusion_dynamics.h"

namespace SPH
{
	/**
	 * @class DiffusionReactionInitialConditionWithRobin
	 * @brief Initialize convection_ for Robin boundary
	 */
	template <class DiffusionReactionParticlesType>
	class DiffusionReactionInitialConditionWithRobin
		: public DiffusionReactionInitialCondition<DiffusionReactionParticlesType>
	{
	public:
		explicit DiffusionReactionInitialConditionWithRobin(SPHBody& sph_body);
		virtual ~DiffusionReactionInitialConditionWithRobin() {};

	protected:
		StdLargeVec<Real> convection_;
		GlobalVariable<Real> T_infinity_;
	};

	/**
	 * @class RelaxationOfAllDiffusionSpeciesRobinContact
	 * @brief Contact diffusion relaxation with Robin boundary condition.
	 */
	template <class ParticlesType, class ContactParticlesType>
	class DiffusionRelaxationRobin
		: public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType>
	{
		StdLargeVec<Vecd>& n_;
		StdVec<StdLargeVec<Real>*> contact_convection_;
		StdVec<Real *> contact_T_infinity_;
		StdVec<StdLargeVec<Vecd>*> contact_n_;

	protected:
		void getDiffusionChangeRateRobinContact(size_t particle_i, size_t particle_j, Real surface_area_ij_Robin, StdLargeVec<Real>& convection_k, Real& T_infinity_k);

	public:
		explicit DiffusionRelaxationRobin(BaseContactRelation& contact_relation);
		virtual ~DiffusionRelaxationRobin() {};

		inline void interaction(size_t index_i, Real dt = 0.0);
	};
}
#endif // PARTICLE_DYNAMICS_DIFFUSION_REACTION_WITH_BOUNDARY_H