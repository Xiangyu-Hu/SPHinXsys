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
#include "interaction_ck.hpp"

namespace SPH
{
template <typename... T>
class UpdateRelation;

template <class ExecutionPolicy, typename... Parameters>
class UpdateRelation<ExecutionPolicy, Inner<Parameters...>>
    : public Interaction<Inner<Parameters...>>, public BaseDynamics<void>
{
  public:
    UpdateRelation(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~UpdateRelation() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    class InteractKernel
        : public Interaction<Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void incrementNeighborSize(UnsignedInt index_i);
        void updateNeighborList(UnsignedInt index_i);

      protected:
        NeighborSearch neighbor_search_;
        Real grid_spacing_squared_;
    };
    typedef UpdateRelation<ExecutionPolicy, Inner<Parameters...>> LocalDynamicsType;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;

    ExecutionPolicy ex_policy_;
    CellLinkedList &cell_linked_list_;
    Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel> kernel_implementation_;
};

template <class ExecutionPolicy, typename... Parameters>
class UpdateRelation<ExecutionPolicy, Contact<Parameters...>>
    : public Interaction<Contact<Parameters...>>, public BaseDynamics<void>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;
    using NeighborCriterion = typename BaseInteraction::InteractKernel::NeighborCriterion;
    using RelationType = Relation<Contact<Parameters...>>;
    using SourceType = typename RelationType::SourceType;
    using TargetType = typename RelationType::TargetType;
    using MaskedCriterion = typename TargetType::template TargetParticleMask<NeighborCriterion>;

  public:
    UpdateRelation(RelationType &contact_relation);
    virtual ~UpdateRelation() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    class InteractKernel
        : public Interaction<Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void incrementNeighborSize(UnsignedInt source_index);
        void updateNeighborList(UnsignedInt source_index);

      protected:
        MaskedCriterion masked_criterion_;
        NeighborSearch neighbor_search_;
    };

    typedef UpdateRelation<ExecutionPolicy, Contact<Parameters...>> LocalDynamicsType;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;
    UniquePtrsKeeper<KernelImplementation> contact_kernel_implementation_ptrs_;
    ExecutionPolicy ex_policy_;
    StdVec<CellLinkedList *> contact_cell_linked_list_;
    StdVec<KernelImplementation *> contact_kernel_implementation_;
};

template <class ExecutionPolicy>
class UpdateRelation<ExecutionPolicy>
{
  public:
    UpdateRelation() {};
    void exec(Real dt = 0.0) {};
};

template <class ExecutionPolicy, class FirstRelation, class... Others>
class UpdateRelation<ExecutionPolicy, FirstRelation, Others...>
    : public UpdateRelation<ExecutionPolicy, FirstRelation>
{
  protected:
    UpdateRelation<ExecutionPolicy, Others...> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit UpdateRelation(
        FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets);
    virtual void exec(Real dt = 0.0) override;
};
} // namespace SPH
#endif // UPDATE_BODY_RELATION_H
