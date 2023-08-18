/**
 * @file 	general_diffusion_reaction_dynamics.hpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_DIFFUSION_REACTION_DYNAMICS_HPP
#define GENERAL_DIFFUSION_REACTION_DYNAMICS_HPP

#include "general_diffusion_reaction_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class ParticlesType>
DiffusionReactionInitialCondition<ParticlesType>::
    DiffusionReactionInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      DiffusionReactionSimpleData<ParticlesType>(sph_body),
      pos_(this->particles_->pos_), all_species_(this->particles_->all_species_) {}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_DIFFUSION_REACTION_DYNAMICS_HPP