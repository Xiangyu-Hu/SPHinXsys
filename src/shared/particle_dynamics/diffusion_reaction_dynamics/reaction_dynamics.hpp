/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	reaction_dynamics.hpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef REACTION_DYNAMICS_HPP
#define REACTION_DYNAMICS_HPP

#include "reaction_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class ParticlesType>
BaseReactionRelaxation<ParticlesType>::
    BaseReactionRelaxation(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      DiffusionReactionSimpleData<ParticlesType>(sph_body),
      reactive_species_(this->particles_->ReactiveSpecies()),
      reaction_model_(this->particles_->diffusion_reaction_material_.ReactionModel()) {}
//=================================================================================================//
template <class ParticlesType>
void BaseReactionRelaxation<ParticlesType>::
    loadLocalSpecies(LocalSpecies &local_species, size_t index_i)
{
    for (size_t k = 0; k != NumReactiveSpecies; ++k)
    {
        local_species[k] = (*reactive_species_[k])[index_i];
    }
}
//=================================================================================================//
template <class ParticlesType>
void BaseReactionRelaxation<ParticlesType>::
    applyGlobalSpecies(LocalSpecies &local_species, size_t index_i)
{
    for (size_t k = 0; k != NumReactiveSpecies; ++k)
    {
        (*reactive_species_[k])[index_i] = local_species[k];
    }
}
//=================================================================================================//
template <class ParticlesType>
Real BaseReactionRelaxation<ParticlesType>::UpdateAReactionSpecies::
operator()(Real input, Real production_rate, Real loss_rate, Real dt) const
{
    return input * exp(-loss_rate * dt) +
           production_rate * (1.0 - exp(-loss_rate * dt)) / (loss_rate + TinyReal);
}
//=================================================================================================//
template <class ParticlesType>
void BaseReactionRelaxation<ParticlesType>::
    advanceForwardStep(size_t index_i, Real dt)
{
    LocalSpecies local_species;
    loadLocalSpecies(local_species, index_i);
    for (size_t k = 0; k != NumReactiveSpecies; ++k)
    {
        Real production_rate = reaction_model_.get_production_rates_[k](local_species);
        Real loss_rate = reaction_model_.get_loss_rates_[k](local_species);
        local_species[k] = updateAReactionSpecies(local_species[k], production_rate, loss_rate, dt);
    }
    applyGlobalSpecies(local_species, index_i);
}
//=================================================================================================//
template <class ParticlesType>
void BaseReactionRelaxation<ParticlesType>::
    advanceBackwardStep(size_t index_i, Real dt)
{
    LocalSpecies local_species;
    loadLocalSpecies(local_species, index_i);
    for (size_t k = NumReactiveSpecies; k != 0; --k)
    {
        size_t m = k - 1;
        Real production_rate = reaction_model_.get_production_rates_[m](local_species);
        Real loss_rate = reaction_model_.get_loss_rates_[m](local_species);
        local_species[m] = updateAReactionSpecies(local_species[m], production_rate, loss_rate, dt);
    }
    applyGlobalSpecies(local_species, index_i);
}
//=================================================================================================//
} // namespace SPH
#endif // REACTION_DYNAMICS_HPP