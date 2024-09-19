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
template <class ReactionModelType>
BaseReactionRelaxation<ReactionModelType>::
    BaseReactionRelaxation(SPHBody &sph_body, ReactionModelType &reaction_model)
    : LocalDynamics(sph_body), reaction_model_(reaction_model)
{
    ReactiveSpeciesNames &species_names = reaction_model.getSpeciesNames();
    for (size_t k = 0; k != NumReactiveSpecies; ++k)
    {
        reactive_species_.push_back(this->particles_->template registerStateVariable<Real>(species_names[k]));
    }
}
//=================================================================================================//
template <class ReactionModelType>
void BaseReactionRelaxation<ReactionModelType>::
    loadLocalSpecies(LocalSpecies &local_species, size_t index_i)
{
    for (size_t k = 0; k != NumReactiveSpecies; ++k)
    {
        local_species[k] = reactive_species_[k][index_i];
    }
}
//=================================================================================================//
template <class ReactionModelType>
void BaseReactionRelaxation<ReactionModelType>::
    applyGlobalSpecies(LocalSpecies &local_species, size_t index_i)
{
    for (size_t k = 0; k != NumReactiveSpecies; ++k)
    {
        reactive_species_[k][index_i] = local_species[k];
    }
}
//=================================================================================================//
template <class ReactionModelType>
Real BaseReactionRelaxation<ReactionModelType>::UpdateReactionSpecies::
operator()(Real input, Real production_rate, Real loss_rate, Real dt) const
{
    Real alpha = exp(-loss_rate * dt);
    return input * alpha + production_rate * (1.0 - alpha) / (loss_rate + TinyReal);
}
//=================================================================================================//
template <class ReactionModelType>
void BaseReactionRelaxation<ReactionModelType>::advanceForwardStep(size_t index_i, Real dt)
{
    LocalSpecies local_species;
    loadLocalSpecies(local_species, index_i);
    for (size_t k = 0; k != NumReactiveSpecies; ++k)
    {
        Real production_rate = reaction_model_.get_production_rates_[k](local_species);
        Real loss_rate = reaction_model_.get_loss_rates_[k](local_species);
        local_species[k] = update_reaction_species_(local_species[k], production_rate, loss_rate, dt);
    }
    applyGlobalSpecies(local_species, index_i);
}
//=================================================================================================//
template <class ReactionModelType>
void BaseReactionRelaxation<ReactionModelType>::advanceBackwardStep(size_t index_i, Real dt)
{
    LocalSpecies local_species;
    loadLocalSpecies(local_species, index_i);
    for (size_t k = NumReactiveSpecies; k != 0; --k)
    {
        size_t m = k - 1;
        Real production_rate = reaction_model_.get_production_rates_[m](local_species);
        Real loss_rate = reaction_model_.get_loss_rates_[m](local_species);
        local_species[m] = update_reaction_species_(local_species[m], production_rate, loss_rate, dt);
    }
    applyGlobalSpecies(local_species, index_i);
}
//=================================================================================================//
} // namespace SPH
#endif // REACTION_DYNAMICS_HPP