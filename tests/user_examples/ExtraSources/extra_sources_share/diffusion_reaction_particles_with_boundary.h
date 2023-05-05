/**
 * @file 	diffusion_reaction_particles_with_boundary.h
 * @brief 	This is the derived class of diffusion reaction particles.
 * @author	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef DIFFUSION_REACTION_PARTICLES_WITH_BOUNDARY_H
#define DIFFUSION_REACTION_PARTICLES_WITH_BOUNDARY_H

#include "base_particles.h"
#include "base_body.h"
#include "base_material.h"
#include "diffusion_reaction.h"
#include "diffusion_reaction_particles.h"

namespace SPH
{
	/**
	 * @class DiffusionReactionParticlesWithBoundary
	 * @brief A group of particles with diffusion or/and reactions particle data.
	 */
	template <class BaseParticlesType, class DiffusionReactionMaterialType>
	class DiffusionReactionParticlesWithBoundary : public DiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType>
	{
	public:
		StdLargeVec<Real> heat_flux_; /**< heat flux for Neumann boundary condition */
		StdLargeVec<Real> convection_; /**< convection for Robin boundary condition */
		StdLargeVec<Real> T_infinity_; /**< T_infinity for Robin boundary condition */

		DiffusionReactionParticlesWithBoundary(SPHBody& sph_body, DiffusionReactionMaterialType* diffusion_reaction_material)
			: DiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType>(sph_body, diffusion_reaction_material)
		{
			heat_flux_.resize(this->all_species_.size());
			convection_.resize(this->all_species_.size());
			T_infinity_.resize(this->all_species_.size());
		};
		virtual ~DiffusionReactionParticlesWithBoundary() {};

		virtual void initializeOtherVariables() override
		{
			DiffusionReactionParticles<BaseParticlesType, DiffusionReactionMaterialType>::initializeOtherVariables();
			
			this->registerVariable(heat_flux_, "HeatFlux");
			this->template addVariableToWrite<Real>("HeatFlux");
			
			this->registerVariable(convection_, "Convection");
			this->template addVariableToWrite<Real>("Convection");

			this->registerVariable(T_infinity_, "T_infinity");
			this->template addVariableToWrite<Real>("T_infinity");
		};

		virtual DiffusionReactionParticlesWithBoundary<BaseParticlesType, DiffusionReactionMaterialType>* ThisObjectPtr() override { return this; };
	};
}

#endif //  DIFFUSION_REACTION_PARTICLES_WITH_BOUNDARY_H