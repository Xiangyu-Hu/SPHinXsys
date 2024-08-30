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
 * @file    interaction_dynamics_ck.h
 * @brief 	This is for the base classes of local particle dynamics, which describe the
 * 			dynamics of a particle and it neighbors.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef INTERACTION_DYNAMICS_CK_H
#define INTERACTION_DYNAMICS_CK_H

#include "base_local_dynamics.h"
#include "neighborhood_ck.hpp"

namespace SPH
{
template <typename... T>
class InteractionDynamics;

template <typename... Parameters>
class InteractionDynamics<Inner<Parameters...>> : public LocalDynamics
{
  public:
    explicit InteractionDynamics(InnerRelation &inner_relation);
    virtual ~InteractionDynamics(){};

    template <class ExecutionPolicy>
    class ComputingKernel : public NeighborList, public Neighbor<Parameters...>
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        InteractionDynamics<Inner<Parameters...>> &encloser);
    };

  protected:
    SPHAdaptation *sph_adaptation_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_index_;
    DiscreteVariable<UnsignedInt> *dv_particle_offset_;
};

template <typename... Parameters>
class InteractionDynamics<Contact<Parameters...>> : public LocalDynamics
{
  public:
    explicit InteractionDynamics(ContactRelation &contact_relation);
    virtual ~InteractionDynamics(){};

    template <class ExecutionPolicy>
    class ComputingKernel : public NeighborList, public Neighbor<Parameters...>
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        InteractionDynamics<Contact<Parameters...>> &encloser,
                        UnsignedInt contact_index);
    };

  protected:
    SPHAdaptation *sph_adaptation_;
    DiscreteVariable<Vecd> *dv_pos_;
    RealBodyVector contact_bodies_;
    StdVec<SPHAdaptation *> contact_adaptations_;
    StdVec<DiscreteVariable<Vecd> *> contact_pos_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_particle_offset_;
};
} // namespace SPH
#endif // INTERACTION_DYNAMICS_CK_H
