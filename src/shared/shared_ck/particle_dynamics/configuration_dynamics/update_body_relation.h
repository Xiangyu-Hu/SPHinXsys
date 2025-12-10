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
#include "base_local_dynamics.h"
#include "base_particles.hpp"
#include "neighborhood_ck.h"
#include "relation_ck.hpp"

namespace SPH
{

template <typename... T>
class UpdateRelation;

template <class ExecutionPolicy, typename... Parameters>
class UpdateRelation<ExecutionPolicy, Inner<Parameters...>>
    : public BaseLocalDynamics<typename Inner<Parameters...>::SourceType>, public BaseDynamics<void>
{
    using BaseLocalDynamicsType = BaseLocalDynamics<typename Inner<Parameters...>::SourceType>;
    using InnerRelationType = Inner<Parameters...>;
    using NeighborList = typename InnerRelationType::NeighborList;
    using Identifier = typename BaseLocalDynamicsType::Identifier;
    using MaskedSource = typename Identifier::SourceParticleMask;
    using NeighborMethodType = typename InnerRelationType::NeighborhoodType;
    using NeighborCriterion = typename NeighborMethodType::NeighborCriterion;
    using MaskedCriterion = typename Identifier::template TargetParticleMask<NeighborCriterion>;

    class OneSidedCheck
    {
        typename NeighborMethodType::ReverseNeighborCriterion reverse_criterion_;

      public:
        template <class EncloserType>
        OneSidedCheck(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : reverse_criterion_(ex_policy, encloser){};
        bool operator()(UnsignedInt i, UnsignedInt j) const
        {
            return i < j || !reverse_criterion_(i, j);
        }
    };

  public:
    UpdateRelation(Inner<Parameters...> &inner_relation);
    virtual ~UpdateRelation() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    class InteractKernel : public NeighborList
    {
      public:
        template <class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void clearNeighborSize(UnsignedInt source_index);
        void incrementNeighborSize(UnsignedInt source_index);
        void updateNeighborList(UnsignedInt source_index);

      protected:
        Vecd *src_pos_;
        UnsignedInt *neighbor_size_;
        MaskedSource masked_src_;
        OneSidedCheck is_one_sided_;
        MaskedCriterion masked_criterion_;
        NeighborSearch neighbor_search_;
    };
    typedef UpdateRelation<ExecutionPolicy, Inner<Parameters...>> LocalDynamicsType;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;

    ExecutionPolicy ex_policy_;
    InnerRelationType &inner_relation_;
    CellLinkedList &cell_linked_list_;
    Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel> kernel_implementation_;
};

template <class ExecutionPolicy, typename... Parameters>
class UpdateRelation<ExecutionPolicy, Contact<Parameters...>>
    : public BaseLocalDynamics<typename Contact<Parameters...>::SourceType>, public BaseDynamics<void>
{
    using BaseLocalDynamicsType = BaseLocalDynamics<typename Contact<Parameters...>::SourceType>;
    using ContactRelationType = Contact<Parameters...>;
    using NeighborList = typename ContactRelationType::NeighborList;
    using Neighborhood = typename ContactRelationType::NeighborhoodType;
    using SearchBox = typename Neighborhood::SearchBox;
    using Identifier = typename BaseLocalDynamicsType::Identifier;
    using SourceType = typename ContactRelationType::SourceType;
    using TargetType = typename ContactRelationType::TargetType;
    using MaskedSource = typename SourceType::SourceParticleMask;
    using NeighborCriterion = typename Neighborhood::NeighborCriterion;
    using MaskedCriterion = typename TargetType::template TargetParticleMask<NeighborCriterion>;

  public:
    UpdateRelation(ContactRelationType &contact_relation);
    virtual ~UpdateRelation() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    class InteractKernel : public NeighborList
    {
      public:
        template <class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void incrementNeighborSize(UnsignedInt source_index);
        void updateNeighborList(UnsignedInt source_index);

      protected:
        Vecd *src_pos_;
        MaskedSource masked_src_;
        MaskedCriterion masked_criterion_;
        NeighborSearch neighbor_search_;
        SearchBox search_box_;
    };

    typedef UpdateRelation<ExecutionPolicy, Contact<Parameters...>> LocalDynamicsType;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;
    UniquePtrsKeeper<KernelImplementation> contact_kernel_implementation_ptrs_;
    ExecutionPolicy ex_policy_;
    ContactRelationType &contact_relation_;
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

template <class ExecutionPolicy, class FirstRelation, class... OtherRelations>
class UpdateRelation<ExecutionPolicy, FirstRelation, OtherRelations...>
    : public UpdateRelation<ExecutionPolicy, FirstRelation>
{
  protected:
    UpdateRelation<ExecutionPolicy, OtherRelations...> other_interactions_;

  public:
    explicit UpdateRelation(FirstRelation &first_relation, OtherRelations &...other_relations);
    virtual void exec(Real dt = 0.0) override;
};
} // namespace SPH
#endif // UPDATE_BODY_RELATION_H
