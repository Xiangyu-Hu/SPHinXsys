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
    virtual ~BodyRelationUpdate(){};

    template <class ExecutionPolicy>
    class ComputingKernel
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        BodyRelationUpdate<Inner<>> &update_inner_relation);
        void incrementNeighborSize(UnsignedInt index_i);
        void updateNeighborList(UnsignedInt index_i);

      protected:
        friend class BodyRelationUpdate<Inner<>>;
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

template <class BodyRelationUpdateType, class ExecutionPolicy>
class UpdateRelation : public BodyRelationUpdateType, public BaseDynamics<void>
{
    using ComputingKernel = typename BodyRelationUpdateType::
        template ComputingKernel<ExecutionPolicy>;

  public:
    template <typename... Args>
    UpdateRelation(Args &&...args);
    virtual ~UpdateRelation(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    ExecutionPolicy ex_policy_;
    Implementation<BodyRelationUpdateType, ExecutionPolicy> kernel_implementation_;
};

template <>
class BodyRelationUpdate<Contact<>> : public LocalDynamics
{

  public:
    explicit BodyRelationUpdate(ContactRelation &contact_relation);
    virtual ~BodyRelationUpdate(){};

    template <class ExecutionPolicy>
    class ComputingKernel
    {
      public:
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        BodyRelationUpdate<Contact<>> &update_contact_relation,
                        UnsignedInt contact_body_index);
        void incrementNeighborSize(UnsignedInt index_i);
        void updateNeighborList(UnsignedInt index_i);

      protected:
        friend class BodyRelationUpdate<Inner<>>;
        NeighborSearch neighbor_search_;
        Vecd *pos_;
        UnsignedInt *neighbor_index_;
        UnsignedInt *particle_offset_;
    };

  protected:
    UnsignedInt particle_offset_list_size_;
    DiscreteVariable<Vecd> *dv_pos_;
    StdVec<CellLinkedList *> contact_cell_linked_list_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_neighbor_index_;
    StdVec<DiscreteVariable<UnsignedInt> *> dv_contact_particle_offset_;
};

template <class BodyRelationUpdateType, class ExecutionPolicy>
class UpdateContactRelation : public BodyRelationUpdateType, public BaseDynamics<void>
{
    using ComputingKernel = typename BodyRelationUpdateType::
        template ComputingKernel<ExecutionPolicy>;
    UniquePtrsKeeper<Implementation<BodyRelationUpdateType, ExecutionPolicy>> contact_kernel_implementation_ptrs_;

  public:
    template <typename... Args>
    UpdateContactRelation(Args &&...args);
    virtual ~UpdateContactRelation(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    ExecutionPolicy ex_policy_;
    StdVec<Implementation<BodyRelationUpdateType, ExecutionPolicy> *> contact_kernel_implementation_;
};
} // namespace SPH
#endif // UPDATE_BODY_RELATION_H
