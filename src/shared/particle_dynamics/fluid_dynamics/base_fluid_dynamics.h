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
 * @file base_fluid_dynamics.h
 * @brief Collection of headers and types used by all fluid dynamics classes.
 * @author Xiangyu Hu
 */

#ifndef BASE_FLUID_DYNAMICS_H
#define BASE_FLUID_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_particles.hpp"
#include "solid_particles.h"

namespace SPH
{
//----------------------------------------------------------------------
// Interaction types specifically for fluid dynamics
//----------------------------------------------------------------------
template <typename InteractionType>
class FreeSurface; /**< A interaction considering the effect of free surface */
template <typename InteractionType>
class FreeStream; /**< A interaction considering the effect of free stream */
template <typename InteractionType>
class AngularConservative; /**< A interaction considering the conservation of angular momentum */

namespace fluid_dynamics
{
typedef DataDelegateSimple<BaseParticles> FluidDataSimple;
typedef DataDelegateInner<BaseParticles> FluidDataInner;
typedef DataDelegateContact<BaseParticles, BaseParticles> FluidContactData;
typedef DataDelegateContact<BaseParticles, SolidParticles, DataDelegateEmptyBase> FluidWallData;
typedef DataDelegateContact<BaseParticles, SolidParticles> FSIContactData;
/**
 * @class InteractionWithWall
 * @brief Base class adding interaction with wall to general relaxation process
 */

template <template <class DataDelegationType> class BaseInteractionType>
class InteractionWithWall : public BaseInteractionType<FSIContactData>
{
  public:
    explicit InteractionWithWall(BaseContactRelation &wall_contact_relation)
        : BaseInteractionType<FSIContactData>(wall_contact_relation)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            wall_vel_ave_.push_back(this->contact_particles_[k]->AverageVelocity());
            wall_acc_ave_.push_back(this->contact_particles_[k]->AverageAcceleration());
            wall_n_.push_back(&(this->contact_particles_[k]->n_));
        }
    };
    virtual ~InteractionWithWall(){};

  protected:
    StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_acc_ave_, wall_n_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // BASE_FLUID_DYNAMICS_H
