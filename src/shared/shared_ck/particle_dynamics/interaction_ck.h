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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file interaction_ck.h
 * @brief This is for the base classes of particle interaction with its neighbors,
 * which describe the basic patterns of particle methods.
 * @author Xiangyu Hu
 */

#ifndef INTERACTION_CK_H
#define INTERACTION_CK_H

#include "base_local_dynamics.h"
#include "neighborhood_ck.h"
#include "relation_ck.hpp"

namespace SPH
{
class WithUpdate;
class WithInitialization;
class OneLevel;
class InteractionOnly;

template <typename... T>
class Interaction;

template <typename... Parameters>
class Interaction<Inner<Parameters...>>
    : public BaseLocalDynamics<typename Inner<Parameters...>::SourceType>
{
    using BaseLocalDynamicsType = BaseLocalDynamics<typename Inner<Parameters...>::SourceType>;
    typedef Inner<Parameters...> InnerRelationType;
    using NeighborList = typename InnerRelationType::NeighborList;
    using Neighborhood = typename InnerRelationType::NeighborhoodType;
    using NeighborKernel = typename Neighborhood::NeighborKernel;

  public:
    explicit Interaction(InnerRelationType &inner_relation);
    virtual ~Interaction() {};

    class InteractKernel : public NeighborList, public NeighborKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
    };

    typedef InteractKernel BaseInteractKernel;
    void registerComputingKernel(Implementation<Base> *implementation);
    void resetComputingKernelUpdated();

  protected:
    InnerRelationType &inner_relation_;
    DiscreteVariable<Real> *dv_Vol_;
};

template <typename... Parameters>
class Interaction<Contact<Parameters...>>
    : public BaseLocalDynamics<typename Contact<Parameters...>::SourceType>
{
    using BaseLocalDynamicsType = BaseLocalDynamics<typename Contact<Parameters...>::SourceType>;
    typedef Contact<Parameters...> ContactRelationType;
    using NeighborList = typename ContactRelationType::NeighborList;
    using Neighborhood = typename ContactRelationType::NeighborhoodType;
    using NeighborKernel = typename Neighborhood::NeighborKernel;

  public:
    explicit Interaction(ContactRelationType &contact_relation);
    virtual ~Interaction() {};

    class InteractKernel : public NeighborList, public NeighborKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                       UnsignedInt contact_index);
    };

    typedef InteractKernel BaseInteractKernel;
    void registerComputingKernel(Implementation<Base> *implementation, UnsignedInt contact_index);
    void resetComputingKernelUpdated(UnsignedInt contact_index);

  protected:
    ContactRelationType &contact_relation_;
    StdVec<SPHBody *> contact_bodies_;
    StdVec<BaseParticles *> contact_particles_;
    StdVec<SPHAdaptation *> contact_adaptations_;
    DiscreteVariable<Real> *dv_Vol_;
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};

template <>
class Interaction<Wall>
{
  public:
    template <class WallContactRelationType>
    Interaction(WallContactRelationType &wall_contact_relation);
    virtual ~Interaction() {};

  protected:
    StdVec<DiscreteVariable<Vecd> *> dv_wall_vel_ave_, dv_wall_acc_ave_, dv_wall_n_;
};
} // namespace SPH
#endif // INTERACTION_CK_H
