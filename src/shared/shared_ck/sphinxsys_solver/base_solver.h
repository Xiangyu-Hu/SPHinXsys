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
 * @file base_solver.h
 * @brief Graph base solver with the help of the Taskflow library.
 * @author	Xiangyu Hu
 */

#ifndef BASE_SOLVER_H
#define BASE_SOLVER_H

#include "base_particle_dynamics.h"
#include "interaction_algorithms_ck.h"
#include "io_base.h"
#include "ownership.h"
#include "particle_sort_ck.h"
#include "simple_algorithms_ck.h"
#include "update_body_relation.h"
#include "update_cell_linked_list.h"

#include <taskflow/taskflow.hpp>

namespace SPH
{
class SPHSolver
{
    UniquePtrsKeeper<BaseDynamics<void>> particle_dynamics_keeper_;
    DataContainerUniquePtrAssemble<BaseDynamics> reduce_dynamics_keeper_;
    UniquePtrKeeper<BodyStatesRecording> state_recording_keeper_;
    UniquePtrsKeeper<BaseIO> io_dynamics_keeper_;

  public:
    SPHSolver(SPHSystem &sph_system);
    virtual ~SPHSolver() {};

    void setRestartIterationStep(size_t iteration_step)
    {
        iteration_step_ = iteration_step;
    };

    template <typename ExecutePolicy, class DynamicsIdentifier, typename... Args>
    auto &addCellLinkedListDynamics(const ExecutePolicy &ex_policy, DynamicsIdentifier &identifier, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            UpdateCellLinkedList<ExecutePolicy, DynamicsIdentifier>>(
            identifier, std::forward<Args>(args)...);
    };

    template <typename ExecutePolicy, class FirstRelation, typename... OtherRelations>
    auto &addRelationDynamics(const ExecutePolicy &ex_policy, FirstRelation &first_relation, OtherRelations &...other_relations)
    {
        return *particle_dynamics_keeper_.createPtr<
            UpdateRelation<ExecutePolicy, FirstRelation, OtherRelations...>>(
            first_relation, other_relations...);
    };

    template <typename ExecutePolicy, typename... Args>
    auto &addSortDynamics(const ExecutePolicy &ex_policy, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            ParticleSortCK<ExecutePolicy>>(std::forward<Args>(args)...);
    };

    template <class UpdateType, typename ExecutePolicy, typename... Args>
    auto &addStateDynamics(const ExecutePolicy &ex_policy, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            StateDynamics<ExecutePolicy, UpdateType>>(std::forward<Args>(args)...);
    };

    template <class ReduceType, typename ExecutePolicy, typename... Args>
    auto &addReduceDynamics(const ExecutePolicy &ex_policy, Args &&...args)
    {
        using FinishDynamics = typename ReduceType::FinishDynamics;
        using OutputType = typename FinishDynamics::OutputType;
        constexpr int type_index = DataTypeIndex<OutputType>::value;
        UniquePtrsKeeper<BaseDynamics<OutputType>> &dynamics_ptrs = std::get<type_index>(reduce_dynamics_keeper_);
        return *dynamics_ptrs.template createPtr<
            ReduceDynamicsCK<ExecutePolicy, ReduceType>>(std::forward<Args>(args)...);
    };

    template <class InteractionType, typename ExecutePolicy, typename... Args>
    auto &addInteractionDynamics(const ExecutePolicy &ex_policy, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<ExecutePolicy, InteractionType>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class InteractionType,
              typename... ControlParameters, typename ExecutePolicy,
              template <typename...> class RelationType, typename... RelationParameters, typename... Args>
    auto &incrementInteractionDynamics(
        const ExecutePolicy &ex_policy, RelationType<RelationParameters...> &inner_relation, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<
                ExecutePolicy, InteractionType<RelationType<ControlParameters..., RelationParameters...>>>>(
            inner_relation, std::forward<Args>(args)...);
    };

    template <class StateRecordingType, typename... Args>
    auto &addStatesRecording(Args &&...args)
    {
        return *state_recording_keeper_.createPtr<
            StateRecordingType>(std::forward<Args>(args)...);
    };

    template <class DynamicsType, typename... Args>
    auto &addIODynamics(Args &&...args)
    {
        return *io_dynamics_keeper_.createPtr<
            DynamicsType>(std::forward<Args>(args)...);
    };

  protected:
    SingularVariable<Real> *physical_time_;
    size_t iteration_step_ = 0;
};
} // namespace SPH
#endif // BASE_SOLVER_H