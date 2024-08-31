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
 * @file    local_interaction_dynamics_ck.h
 * @brief 	This is for the base classes of local particle dynamics, which describe the
 * 			dynamics of a particle and it neighbors.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef LOCAL_INTERACTION_DYNAMICS_CK_H
#define LOCAL_INTERACTION_DYNAMICS_CK_H

#include "base_local_dynamics.h"
#include "neighborhood_ck.hpp"

namespace SPH
{
template <typename... T>
class LocalInteractionDynamics;

template <typename... Parameters>
class LocalInteractionDynamics<Inner<Parameters...>> : public LocalDynamics
{
  public:
    explicit LocalInteractionDynamics(InnerRelation &inner_relation);
    virtual ~LocalInteractionDynamics(){};

    class ComputingKernel : public NeighborList, public Neighbor<Parameters...>
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        LocalInteractionDynamics<Inner<Parameters...>> &encloser);
    };

    void registerComputingKernel(Implementation<Base> *implementation);
    void resetComputingKernelUpdated();

  protected:
    InnerRelation &inner_relation_;
    SPHAdaptation *sph_adaptation_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_index_;
    DiscreteVariable<UnsignedInt> *dv_particle_offset_;
};

template <typename... Parameters>
class LocalInteractionDynamics<Contact<Parameters...>> : public LocalDynamics
{
  public:
    explicit LocalInteractionDynamics(ContactRelation &contact_relation);
    virtual ~LocalInteractionDynamics(){};

    class ComputingKernel : public NeighborList, public Neighbor<Parameters...>
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        LocalInteractionDynamics<Contact<Parameters...>> &encloser,
                        UnsignedInt contact_index);
    };

    void registerComputingKernel(Implementation<Base> *implementation, UnsignedInt contact_index);
    void resetComputingKernelUpdated(UnsignedInt contact_index);

  protected:
    ContactRelation &contact_relation_;
    SPHAdaptation *sph_adaptation_;
    DiscreteVariable<Vecd> *dv_pos_;
    RealBodyVector contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    StdVec<SPHAdaptation *> contact_adaptations_;
    StdVec<DiscreteVariable<Vecd> *> contact_pos_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_particle_offset_;
};
} // namespace SPH
#endif // LOCAL_INTERACTION_DYNAMICS_CK_H
