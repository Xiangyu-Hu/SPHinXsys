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
 * @file particle_method_container.h
 * @brief interface for simplify the modeling process,
 * especially decrease the explicit usage of template by type deduction.
 * @author	Xiangyu Hu
 */

#ifndef PARTICLE_METHOD_CONTAINER_H
#define PARTICLE_METHOD_CONTAINER_H

#include "base_particle_dynamics.h"
#include "interaction_algorithms_ck.h"
#include "io_base.h"
#include "ownership.h"
#include "particle_sort_ck.h"
#include "simple_algorithms_ck.h"
#include "update_body_relation.h"
#include "update_cell_linked_list.h"

namespace SPH
{

class BaseMethodContainer
{
  public:
    virtual ~BaseMethodContainer(){};
};

template <typename ExecutePolicy>
class ParticleMethodContainer : public BaseMethodContainer
{
    UniquePtrsKeeper<BaseDynamics<void>> particle_dynamics_keeper_;
    DataContainerUniquePtrAssemble<BaseDynamics> reduce_dynamics_keeper_;
    UniquePtrsKeeper<BodyStatesRecording> state_recorders_keeper_;
    UniquePtrsKeeper<BaseIO> other_io_keeper_;

  public:
    ParticleMethodContainer(const ExecutePolicy &ex_policy)
    : BaseMethodContainer() {};
    virtual ~ParticleMethodContainer() {};

    template <class DynamicsIdentifier, typename... Args>
    auto &addCellLinkedListDynamics(DynamicsIdentifier &identifier, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            UpdateCellLinkedList<ExecutePolicy, DynamicsIdentifier>>(
            identifier, std::forward<Args>(args)...);
    };

    template <class FirstRelation, typename... OtherRelations>
    auto &addRelationDynamics(FirstRelation &first_relation, OtherRelations &...other_relations)
    {
        return *particle_dynamics_keeper_.createPtr<
            UpdateRelation<ExecutePolicy, FirstRelation, OtherRelations...>>(
            first_relation, other_relations...);
    };

    template <typename... Args>
    auto &addSortDynamics(Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            ParticleSortCK<ExecutePolicy>>(std::forward<Args>(args)...);
    };

    template <class UpdateType, typename... Args>
    auto &addStateDynamics(Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            StateDynamics<ExecutePolicy, UpdateType>>(std::forward<Args>(args)...);
    };

    template <class ReduceType, typename... Args>
    auto &addReduceDynamics(Args &&...args)
    {
        using FinishDynamics = typename ReduceType::FinishDynamics;
        using OutputType = typename FinishDynamics::OutputType;
        constexpr int type_index = DataTypeIndex<OutputType>::value;
        UniquePtrsKeeper<BaseDynamics<OutputType>> &dynamics_ptrs = std::get<type_index>(reduce_dynamics_keeper_);
        return *dynamics_ptrs.template createPtr<
            ReduceDynamicsCK<ExecutePolicy, ReduceType>>(std::forward<Args>(args)...);
    };

    template <class InteractionType, typename... Args>
    auto &addInteractionDynamics(Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<ExecutePolicy, InteractionType>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class InteractionType,
              typename... ControlParameters,
              template <typename...> class RelationType, typename... RelationParameters, typename... Args>
    auto &addInteractionDynamics(
        RelationType<RelationParameters...> &relation, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<
                ExecutePolicy, InteractionType<RelationType<ControlParameters..., RelationParameters...>>>>(
            relation, std::forward<Args>(args)...);
    };

    template <template <typename...> class RecorderType, typename... Args>
    auto &addBodyStateRecorder(Args &&...args)
    {
        return *state_recorders_keeper_.createPtr<RecorderType<ExecutePolicy>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class RegressionType,
              template <typename...> class ObservationType, typename... Parameters, typename... Args>
    auto &addRegressionTest(Args &&...args)
    {
        return *other_io_keeper_.createPtr<
            RegressionType<ObservationType<ExecutePolicy, Parameters...>>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class IOType, typename... Parameters, typename... Args>
    auto &addIODynamics(Args &&...args)
    {
        return *other_io_keeper_.createPtr<
            IOType<ExecutePolicy, Parameters...>>(std::forward<Args>(args)...);
    };
};
} // namespace SPH
#endif // PARTICLE_METHOD_CONTAINER_H