/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file reaction_dynamics.h
 * @brief Opertor spliting method for unconditionally stable time stepping.
 * @author Chi Zhang and Xiangyu Hu
 */

#ifndef REACTION_DYNAMICS_H
#define REACTION_DYNAMICS_H

#include "general_diffusion_reaction_dynamics.h"

namespace SPH
{
/**
 * @class BaseReactionRelaxation
 * @brief Base class for computing the reaction process of all species
 */
template <class ReactionModelType>
class BaseReactionRelaxation : public LocalDynamics
{
  protected:
    struct UpdateReactionSpecies
    {
        Real operator()(Real input, Real production_rate, Real loss_rate, Real dt) const;
    };
    void advanceForwardStep(size_t index_i, Real dt);
    void advanceBackwardStep(size_t index_i, Real dt);

  private:
    static constexpr int NumReactiveSpecies = ReactionModelType::NumSpecies;
    typedef std::array<std::string, NumReactiveSpecies> ReactiveSpeciesNames;
    typedef std::array<Real, NumReactiveSpecies> LocalSpecies;
    StdVec<Real *> reactive_species_;
    ReactionModelType &reaction_model_;
    UpdateReactionSpecies update_reaction_species_;
    void loadLocalSpecies(LocalSpecies &local_species, size_t index_i);
    void applyGlobalSpecies(LocalSpecies &local_species, size_t index_i);

  public:
    explicit BaseReactionRelaxation(SPHBody &sph_body, ReactionModelType &reaction_model);
    virtual ~BaseReactionRelaxation(){};
};

/**
 * @class ReactionRelaxationForward
 * @brief Compute the reaction process of all species by forward splitting
 */
template <class ReactionModelType>
class ReactionRelaxationForward : public BaseReactionRelaxation<ReactionModelType>
{
  public:
    template <typename... Args>
    ReactionRelaxationForward(Args &&...args)
        : BaseReactionRelaxation<ReactionModelType>(std::forward<Args>(args)...){};
    virtual ~ReactionRelaxationForward(){};
    void update(size_t index_i, Real dt = 0.0) { this->advanceForwardStep(index_i, dt); };
};

/**
 * @class ReactionRelaxationBackward
 * @brief Compute the reaction process of all species by backward splitting
 */
template <class ReactionModelType>
class ReactionRelaxationBackward : public BaseReactionRelaxation<ReactionModelType>
{
  public:
    template <typename... Args>
    ReactionRelaxationBackward(Args &&...args)
        : BaseReactionRelaxation<ReactionModelType>(std::forward<Args>(args)...){};
    virtual ~ReactionRelaxationBackward(){};
    void update(size_t index_i, Real dt = 0.0) { this->advanceBackwardStep(index_i, dt); };
};
} // namespace SPH
#endif // REACTION_DYNAMICS_H