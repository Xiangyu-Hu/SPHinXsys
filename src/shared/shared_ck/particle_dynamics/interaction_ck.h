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
 * @file    interaction_ck.h
 * @brief 	This is for the base classes of local particle dynamics, which describe the
 * 			dynamics of a particle and it neighbors.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef INTERACTION_CK_H
#define INTERACTION_CK_H

#include "base_local_dynamics.h"
#include "neighborhood_ck.hpp"
#include "relation_ck.hpp"

namespace SPH
{
class WithUpdate;
class WithInitialization;
class OneLevel;

template <typename... T>
class Interaction;

template <typename... Parameters>
class Interaction<Inner<Parameters...>> : public LocalDynamics
{
  public:
    explicit Interaction(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~Interaction(){};

    class InteractKernel : public NeighborList, public Neighbor<Parameters...>
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       Interaction<Inner<Parameters...>> &encloser);
    };

    void registerComputingKernel(Implementation<Base> *implementation);
    void resetComputingKernelUpdated();

  protected:
    Relation<Inner<Parameters...>> &inner_relation_;
    SPHAdaptation *sph_adaptation_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_index_;
    DiscreteVariable<UnsignedInt> *dv_particle_offset_;
};

template <typename... Parameters>
class Interaction<Contact<Parameters...>> : public LocalDynamics
{
  public:
    explicit Interaction(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~Interaction(){};

    class InteractKernel : public NeighborList, public Neighbor<Parameters...>
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       Interaction<Contact<Parameters...>> &encloser,
                       UnsignedInt contact_index);
    };

    void registerComputingKernel(Implementation<Base> *implementation, UnsignedInt contact_index);
    void resetComputingKernelUpdated(UnsignedInt contact_index);

  protected:
    Relation<Contact<Parameters...>> &contact_relation_;
    SPHAdaptation *sph_adaptation_;
    DiscreteVariable<Vecd> *dv_pos_;
    RealBodyVector contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    StdVec<SPHAdaptation *> contact_adaptations_;
    StdVec<DiscreteVariable<Vecd> *> contact_pos_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_particle_offset_;
};

template <typename... Parameters>
class Interaction<Contact<Wall, Parameters...>> : public Interaction<Contact<Parameters...>>
{
  public:
    explicit Interaction(Relation<Contact<Parameters...>> &wall_contact_relation);
    virtual ~Interaction(){};

  protected:
    StdVec<DiscreteVariable<Vecd> *> dv_wall_vel_ave_, dv_wall_acc_ave_, dv_wall_n_;
    StdVec<DiscreteVariable<Real> *> dv_wall_Vol_;
};
} // namespace SPH
#endif // INTERACTION_CK_H
