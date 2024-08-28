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
 * @file    update_body_relation.h
 * @brief   Collection of dynamics for particle configuration.
 * @author	Xiangyu Hu
 */

#ifndef UPDATE_BODY_RELATION_H
#define UPDATE_BODY_RELATION_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_configuration_dynamics.h"
#include "base_particles.hpp"

namespace SPH
{
template <typename... T>
class BodyRelationUpdate;

template <>
class BodyRelationUpdate<Inner<>> : public LocalDynamics
{

  public:
    explicit BodyRelationUpdate(InnerRelation &inner_relation);
    virtual ~BodyRelationUpdate() {};

    template <class ExecutionPolicy>
    class ComputingKernel
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        BodyRelationUpdate<Inner<>> &encloser);
        void incrementNeighborSize(UnsignedInt index_i);
        void updateNeighborList(UnsignedInt index_i);

      protected:
        NeighborSearch neighbor_search_;
        Vecd *pos_;
        UnsignedInt *neighbor_index_;
        UnsignedInt *particle_offset_;
    };

  protected:
    CellLinkedList &cell_linked_list_;
    UnsignedInt particle_offset_list_size_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<UnsignedInt> *dv_neighbor_index_;
    DiscreteVariable<UnsignedInt> *dv_particle_offset_;
};

template <typename... T>
class UpdateRelation;

template <class ExecutionPolicy, typename... Parameters>
class UpdateRelation<ExecutionPolicy, BodyRelationUpdate<Inner<Parameters...>>>
    : public BodyRelationUpdate<Inner<Parameters...>>, public BaseDynamics<void>
{
    typedef BodyRelationUpdate<Inner<Parameters...>> LocalDynamicsType;
    using ComputingKernel = typename LocalDynamicsType::template ComputingKernel<ExecutionPolicy>;
    using KernelImplementation = Implementation<LocalDynamicsType, ExecutionPolicy>;

  public:
    template <typename... Args>
    UpdateRelation(Args &&...args);
    virtual ~UpdateRelation() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    ExecutionPolicy ex_policy_;
    KernelImplementation kernel_implementation_;
};

template <>
class BodyRelationUpdate<Contact<>> : public LocalDynamics
{

  public:
    explicit BodyRelationUpdate(ContactRelation &contact_relation);
    virtual ~BodyRelationUpdate() {};

    template <class ExecutionPolicy>
    class ComputingKernel
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        BodyRelationUpdate<Contact<>> &encloser,
                        UnsignedInt contact_index);
        void incrementNeighborSize(UnsignedInt index_i);
        void updateNeighborList(UnsignedInt index_i);

      protected:
        NeighborSearch neighbor_search_;
        Vecd *pos_;
        UnsignedInt *neighbor_index_;
        UnsignedInt *particle_offset_;
    };

  protected:
    RealBodyVector contact_bodies_;
    UnsignedInt particle_offset_list_size_;
    DiscreteVariable<Vecd> *dv_pos_;
    StdVec<CellLinkedList *> contact_cell_linked_list_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_particle_offset_;
};

template <class ExecutionPolicy, typename... Parameters>
class UpdateRelation<ExecutionPolicy, BodyRelationUpdate<Contact<Parameters...>>>
    : public BodyRelationUpdate<Contact<Parameters...>>, public BaseDynamics<void>
{
    typedef BodyRelationUpdate<Contact<Parameters...>> LocalDynamicsType;
    using ComputingKernel = typename LocalDynamicsType::template ComputingKernel<ExecutionPolicy>;
    using KernelImplementation = Implementation<LocalDynamicsType, ExecutionPolicy>;
    UniquePtrsKeeper<KernelImplementation> contact_kernel_implementation_ptrs_;

  public:
    template <typename... Args>
    UpdateRelation(Args &&...args);
    virtual ~UpdateRelation() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    ExecutionPolicy ex_policy_;
    StdVec<KernelImplementation *> contact_kernel_implementation_;
};

} // namespace SPH
#endif // UPDATE_BODY_RELATION_H
