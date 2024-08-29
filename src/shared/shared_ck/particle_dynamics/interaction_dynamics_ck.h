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
#include "neighborhood.hpp"

namespace SPH
{
template <typename... T>
class InteractionDynamics;

template <>
class InteractionDynamics<Inner<>> : public LocalDynamics
{
  public:
    explicit InteractionDynamics(InnerRelation &inner_relation);
    virtual ~InteractionDynamics(){};

    template <class ExecutionPolicy>
    class ComputingKernel : public NeighborList
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        InteractionDynamics<Inner<>> &encloser);

      protected:
        Vecd *pos_;
        SmoothingKernel kernel_;

        Vecd r_ij(size_t index_i, size_t index_j) const { return pos_[index_i] - pos_[index_j]; };
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    SmoothingKernel *kernel_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_index_;
    DiscreteVariable<UnsignedInt> *dv_particle_offset_;
};

template <>
class InteractionDynamics<Contact<>> : public LocalDynamics
{
  public:
    explicit InteractionDynamics(InnerRelation &inner_relation);
    virtual ~InteractionDynamics(){};

    template <class ExecutionPolicy>
    class ComputingKernel : public NeighborList
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        InteractionDynamics<Inner<>> &encloser,
                        UnsignedInt contact_index);

      protected:
        Vecd *pos_;
        Vecd *contact_pos_;
        SmoothingKernel kernel_;

        Vecd r_ij(size_t index_i, size_t index_j) const { return pos_[index_i] - contact_pos_[index_j]; };
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    StdVec<DiscreteVariable<Vecd> *> contact_pos_;
    StdVec<SmoothingKernel *> contact_kernel_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_particle_offset_;
};
} // namespace SPH
#endif // INTERACTION_DYNAMICS_CK_H
