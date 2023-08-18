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
 * @file    general_diffusion_reaction_dynamics.h
 * @brief   This is the particle dynamics applicable for all type of particles.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_DIFFUSION_REACTION_DYNAMICS_H
#define GENERAL_DIFFUSION_REACTION_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "diffusion_reaction.h"
#include "diffusion_reaction_particles.h"

namespace SPH
{
template <class ParticlesType>
using DiffusionReactionSimpleData = DataDelegateSimple<ParticlesType>;

template <class ParticlesType>
using DiffusionReactionInnerData = DataDelegateInner<ParticlesType>;

template <class ParticlesType, class ContactParticlesType>
using DiffusionReactionContactData =
    DataDelegateContact<ParticlesType, ContactParticlesType>;

/**
 * @class DiffusionReactionInitialCondition
 * @brief Pure abstract class for initial conditions
 */
template <class ParticlesType>
class DiffusionReactionInitialCondition
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    explicit DiffusionReactionInitialCondition(SPHBody &sph_body);
    virtual ~DiffusionReactionInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdVec<StdLargeVec<Real>> &all_species_;
};

/**
 * @class DiffusionReactionSpeciesConstraint
 * @brief set boundary condition for diffusion problem
 */
template <class DynamicsIdentifier, class ParticlesType>
class DiffusionReactionSpeciesConstraint
    : public BaseLocalDynamics<DynamicsIdentifier>,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    DiffusionReactionSpeciesConstraint(DynamicsIdentifier &identifier, const std::string &species_name)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          DiffusionReactionSimpleData<ParticlesType>(identifier.getSPHBody()),
          phi_(this->particles_->diffusion_reaction_material_.AllSpeciesIndexMap()[species_name]),
          species_(this->particles_->all_species_[phi_]){};
    virtual ~DiffusionReactionSpeciesConstraint(){};

  protected:
    size_t phi_;
    StdLargeVec<Real> &species_;
};

/**
 * @class DiffusionBasedMapping
 * @brief Mapping inside of body according to diffusion.
 * This is a abstract class to be override for case specific implementation
 */
template <class ParticlesType>
class DiffusionBasedMapping
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    explicit DiffusionBasedMapping(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          DiffusionReactionSimpleData<ParticlesType>(sph_body),
          pos_(this->particles_->pos_), all_species_(this->particles_->all_species_){};
    virtual ~DiffusionBasedMapping(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdVec<StdLargeVec<Real>> &all_species_;
};

/**
 * @class   SpeciesSummation
 * @brief   Computing the total averaged parameter on the whole diffusion body.
 */
template <class DynamicsIdentifier, class ParticlesType>
class SpeciesSummation
    : public BaseLocalDynamicsReduce<Real, ReduceSum<Real>, DynamicsIdentifier>,
      public DiffusionReactionSimpleData<ParticlesType>
{
  protected:
    StdVec<StdLargeVec<Real>> &all_species_;
    size_t phi_;

  public:
    SpeciesSummation(DynamicsIdentifier &identifier, const std::string &species_name)
        : BaseLocalDynamicsReduce<Real, ReduceSum<Real>, DynamicsIdentifier>(identifier, Real(0)),
          DiffusionReactionSimpleData<ParticlesType>(identifier.getSPHBody()),
          all_species_(this->particles_->all_species_),
          phi_(this->particles_->diffusion_reaction_material_.AllSpeciesIndexMap()[species_name])
    {
        this->quantity_name_ = "DiffusionReactionSpeciesAverage";
    };
    virtual ~SpeciesSummation(){};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return all_species_[phi_][index_i];
    };
};
} // namespace SPH
#endif // GENERAL_DIFFUSION_REACTION_DYNAMICS_H