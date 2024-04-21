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
class FreeSurface;         /**< A interaction considering the effect of free surface */
class FreeStream;          /**< A interaction considering the effect of free stream */
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

template <template <typename...> class BaseInteractionType>
class InteractionWithWall : public BaseInteractionType<FSIContactData>
{
  public:
    explicit InteractionWithWall(BaseContactRelation &wall_contact_relation)
        : BaseInteractionType<FSIContactData>(wall_contact_relation),
          wall_mass_device_(this->contact_particles_.size(), executionQueue.getQueue()),
          wall_Vol_device_(this->contact_particles_.size(), executionQueue.getQueue()),
          wall_vel_ave_device_(this->contact_particles_.size(), executionQueue.getQueue()),
          wall_force_ave_device_(this->contact_particles_.size(), executionQueue.getQueue()),
          wall_n_device_(this->contact_particles_.size(), executionQueue.getQueue())
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            wall_vel_ave_.push_back(this->contact_particles_[k]->AverageVelocity());
            wall_force_ave_.push_back(this->contact_particles_[k]->AverageForce());
            wall_n_.push_back(&(this->contact_particles_[k]->n_));
            wall_mass_.push_back(&(this->contact_particles_[k]->mass_));
            wall_Vol_.push_back(&(this->contact_particles_[k]->Vol_));

            // Device variables
            wall_mass_device_.at(k) = this->contact_particles_[k]->template getDeviceVariableByName<DeviceReal>("Mass");;
            wall_Vol_device_.at(k) = this->contact_particles_[k]->template getDeviceVariableByName<DeviceReal>("Volume");;
            wall_vel_ave_device_.at(k) = this->contact_particles_[k]->template getDeviceVariableByName<DeviceVecd>("Velocity");
            wall_force_ave_device_.at(k) = this->contact_particles_[k]->template getDeviceVariableByName<DeviceVecd>("Force");
            wall_n_device_.at(k) = this->contact_particles_[k]->template getDeviceVariableByName<DeviceVecd>("Normal");
        }
    };
    virtual ~InteractionWithWall(){};

  protected:
    StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_force_ave_, wall_n_;
    StdVec<StdLargeVec<Real> *> wall_mass_, wall_Vol_;

    StdSharedVec<DeviceReal*> wall_mass_device_, wall_Vol_device_;
    StdSharedVec<DeviceVecd*> wall_vel_ave_device_, wall_force_ave_device_, wall_n_device_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // BASE_FLUID_DYNAMICS_H
