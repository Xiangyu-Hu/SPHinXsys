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
 * @file particle_method_container.h
 * @brief interface for simplify the modeling process,
 * especially decrease the explicit usage of template by type deduction.
 * @author	Xiangyu Hu
 */

#ifndef PARTICLE_METHOD_CONTAINER_H
#define PARTICLE_METHOD_CONTAINER_H

#include "base_particle_dynamics.h"
#include "complex_algorithms_ck.h"
#include "interaction_algorithms_ck.h"
#include "io_base.h"
#include "io_observation_ck.h"
#include "ownership.h"
#include "particle_sort_ck.h"
#include "simple_algorithms_ck.h"
#include "update_body_relation.h"
#include "update_cell_linked_list.h"

namespace SPH
{

class ParticleDynamicsGroup : public BaseDynamics<void>
{
    StdVec<BaseDynamics<void> *> particle_dynamics_;

  public:
    ParticleDynamicsGroup() : BaseDynamics<void>() {};
    ParticleDynamicsGroup(const StdVec<BaseDynamics<void> *> &particle_dynamics)
        : BaseDynamics<void>(), particle_dynamics_(particle_dynamics) {}
    ~ParticleDynamicsGroup() {};

    void add(BaseDynamics<void> *dynamics)
    {
        particle_dynamics_.push_back(dynamics);
    }

    void add(const ParticleDynamicsGroup &dynamics)
    {
        for (auto *dynamics : dynamics.getAllDynamics())
        {
            particle_dynamics_.push_back(dynamics);
        }
    }

    ParticleDynamicsGroup operator+(const ParticleDynamicsGroup &other) const
    {
        StdVec<BaseDynamics<void> *> other_dynamics = other.getAllDynamics();
        StdVec<BaseDynamics<void> *> result_dynamics = particle_dynamics_;
        result_dynamics.insert(result_dynamics.end(), other_dynamics.begin(), other_dynamics.end());
        return ParticleDynamicsGroup(result_dynamics);
    }

    StdVec<BaseDynamics<void> *> getAllDynamics() const
    {
        return particle_dynamics_;
    }

    void exec(Real dt = 0.0) override
    {
        for (UnsignedInt i = 0; i != particle_dynamics_.size(); ++i)
        {
            particle_dynamics_[i]->exec(dt);
        }
    }
};

template <class Operation>
class ReduceDynamicsGroup : public BaseDynamics<typename Operation::ReturnType>
{
    Operation operation_;
    using ReturnType = typename Operation::ReturnType;
    StdVec<BaseDynamics<ReturnType> *> reduce_dynamics_;

  public:
    ReduceDynamicsGroup() : BaseDynamics<ReturnType>(), operation_() {}
    ReduceDynamicsGroup(const Operation &operation)
        : BaseDynamics<ReturnType>(), operation_(operation) {}
    ReduceDynamicsGroup(const Operation &operation,
                        const StdVec<BaseDynamics<ReturnType> *> &reduce_dynamics)
        : BaseDynamics<ReturnType>(),
          operation_(operation), reduce_dynamics_(reduce_dynamics) {}
    ~ReduceDynamicsGroup() = default;

    void add(BaseDynamics<ReturnType> *dynamics)
    {
        reduce_dynamics_.push_back(dynamics);
    }

    void add(const ReduceDynamicsGroup<Operation> &dynamics)
    {
        for (auto *dynamics : dynamics.getAllDynamics())
        {
            reduce_dynamics_.push_back(dynamics);
        }
    }

    template <class DerivedReduceDynamicsType>
    void add(const StdVec<DerivedReduceDynamicsType *> &dynamics)
    {
        for (auto *dynamics : dynamics)
        {
            reduce_dynamics_.push_back(dynamics);
        }
    }

    StdVec<BaseDynamics<ReturnType> *> getAllDynamics() const
    {
        return reduce_dynamics_;
    }

    ReduceDynamicsGroup<Operation> operator+(
        const ReduceDynamicsGroup<Operation> &other) const
    {
        StdVec<BaseDynamics<ReturnType> *> other_dynamics = other.getAllDynamics();
        StdVec<BaseDynamics<ReturnType> *> result_dynamics = reduce_dynamics_;
        result_dynamics.insert(result_dynamics.end(), other_dynamics.begin(), other_dynamics.end());
        return ReduceDynamicsGroup<Operation>(operation_, result_dynamics);
    }

    ReturnType exec(Real dt = 0.0) override
    {
        ReturnType result = ReduceReference<Operation>::value;
        for (auto *dynamics : reduce_dynamics_)
        {
            result = operation_(result, dynamics->exec(dt));
        }
        return result;
    };
};

class BaseMethodContainer
{
  public:
    virtual ~BaseMethodContainer() {};
};

template <typename ExecutionPolicy>
class ParticleMethodContainer : public BaseMethodContainer
{
    UniquePtrsKeeper<AbstractDynamics> particle_dynamics_keeper_;
    UniquePtrsKeeper<BodyStatesRecording> state_recorders_keeper_;
    UniquePtrsKeeper<BaseIO> other_io_keeper_;

  public:
    ParticleMethodContainer(const ExecutionPolicy &ex_policy)
        : BaseMethodContainer() {};
    virtual ~ParticleMethodContainer() {};

    template <template <typename...> class GeneralDynamicsType, typename... Parameters, class DynamicsIdentifier, typename... Args>
    auto &addGeneralDynamics(DynamicsIdentifier &identifier, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            GeneralDynamicsType<ExecutionPolicy, Parameters...>>(identifier, std::forward<Args>(args)...);
    };

    template <class GeneralDynamicsType, class DynamicsIdentifier, typename... Args>
    auto &addReturnDynamics(DynamicsIdentifier &identifier, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<GeneralDynamicsType>(identifier, std::forward<Args>(args)...);
    };

    template <class GeneralDynamicsType, class DynamicsIdentifier, typename... Args>
    StdVec<GeneralDynamicsType *> addReturnDynamics(StdVec<DynamicsIdentifier *> &identifiers, Args &&...args)
    {
        StdVec<GeneralDynamicsType *> result;
        for (auto *identifier : identifiers)
        {
            result.push_back(&addReturnDynamics<GeneralDynamicsType>(*identifier, std::forward<Args>(args)...));
        }
        return result;
    };

    template <class DynamicsIdentifier, typename... Args>
    auto &addCellLinkedListDynamics(DynamicsIdentifier &identifier, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            UpdateCellLinkedList<ExecutionPolicy, DynamicsIdentifier>>(
            identifier, std::forward<Args>(args)...);
    };

    template <class DynamicsIdentifier>
    ParticleDynamicsGroup addCellLinkedListDynamics(StdVec<DynamicsIdentifier *> &identifiers)
    {
        ParticleDynamicsGroup group;
        for (auto *identifier : identifiers)
        {
            group.add(&addCellLinkedListDynamics(*identifier));
        }
        return group;
    }

    template <class FirstRelation, typename... OtherRelations>
    auto &addRelationDynamics(FirstRelation &first_relation, OtherRelations &...other_relations)
    {
        return *particle_dynamics_keeper_.createPtr<
            UpdateRelation<ExecutionPolicy, FirstRelation, OtherRelations...>>(
            first_relation, other_relations...);
    };

    template <typename... Args>
    auto &addSortDynamics(Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            ParticleSortCK<ExecutionPolicy>>(std::forward<Args>(args)...);
    };

    template <class DynamicsIdentifier>
    ParticleDynamicsGroup addSortDynamics(StdVec<DynamicsIdentifier *> &identifiers)
    {
        ParticleDynamicsGroup group;
        for (auto *identifier : identifiers)
        {
            group.add(&addSortDynamics(*identifier));
        }
        return group;
    }

    template <class UpdateType, class DynamicsIdentifier, typename... Args>
    auto &addStateDynamics(DynamicsIdentifier &dynamics_identifier, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            StateDynamics<ExecutionPolicy, UpdateType>>(dynamics_identifier, std::forward<Args>(args)...);
    };

    template <class UpdateType, class DynamicsIdentifier, typename... Args>
    ParticleDynamicsGroup addStateDynamics(StdVec<DynamicsIdentifier *> &identifiers, Args &&...args)
    {
        ParticleDynamicsGroup group;
        for (auto &identifier : identifiers)
        {
            group.add(&addStateDynamics<UpdateType>(*identifier, std::forward<Args>(args)...));
        }
        return group;
    }

    template <template <typename...> class UpdateType, typename... ControlParameters,
              class DynamicsIdentifier, typename... Args>
    auto &addStateDynamics(DynamicsIdentifier &dynamics_identifier, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            StateDynamics<ExecutionPolicy, UpdateType<ControlParameters..., DynamicsIdentifier>>>(
            dynamics_identifier, std::forward<Args>(args)...);
    };

    template <template <typename...> class UpdateType, typename... ControlParameters,
              class DynamicsIdentifier, typename... Args>
    ParticleDynamicsGroup addStateDynamics(StdVec<DynamicsIdentifier *> &identifiers, Args &&...args)
    {
        ParticleDynamicsGroup group;
        for (auto &identifier : identifiers)
        {
            group.add(&addStateDynamics<UpdateType, ControlParameters...>(*identifier, std::forward<Args>(args)...));
        }
        return group;
    }

    template <class ReduceType, typename... Args>
    auto &addReduceDynamics(Args &&...args)
    {
        return *particle_dynamics_keeper_.template createPtr<
            ReduceDynamicsCK<ExecutionPolicy, ReduceType>>(std::forward<Args>(args)...);
    };

    template <typename Operation, class ReduceType, typename DynamicsIdentifier, typename... Args>
    ReduceDynamicsGroup<Operation> addReduceDynamics(const StdVec<DynamicsIdentifier *> &identifiers, Args &&...args)
    {
        StdVec<BaseDynamics<typename Operation::ReturnType> *> reduce_dynamics;
        for (auto &identifier : identifiers)
        {
            reduce_dynamics.push_back(&addReduceDynamics<ReduceType>(*identifier, std::forward<Args>(args)...));
        }
        return ReduceDynamicsGroup<Operation>(Operation(), reduce_dynamics);
    };

    template <class InteractionType, typename... Args>
    auto &addInteractionDynamics(Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<ExecutionPolicy, InteractionType>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class InteractionType, typename... ControlParameters,
              template <typename...> class RelationType, typename... RelationParameters, typename... Args>
    auto &addInteractionDynamics(
        RelationType<RelationParameters...> &relation, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<
                ExecutionPolicy, InteractionType<RelationType<ControlParameters..., RelationParameters...>>>>(
            relation, std::forward<Args>(args)...);
    };

    template <template <typename...> class InteractionType, typename... ControlParameters,
              template <typename...> class RelationType, typename... RelationParameters, typename... Args>
    auto &addInteractionDynamicsOneLevel(
        RelationType<RelationParameters...> &relation, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<
                ExecutionPolicy, InteractionType<RelationType<OneLevel, ControlParameters..., RelationParameters...>>>>(
            relation, std::forward<Args>(args)...);
    };

    template <template <typename...> class InteractionType, typename... ControlParameters,
              template <typename...> class RelationType, typename... RelationParameters, typename... Args>
    auto &addInteractionDynamicsWithUpdate(
        RelationType<RelationParameters...> &relation, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<
                ExecutionPolicy, InteractionType<RelationType<WithUpdate, ControlParameters..., RelationParameters...>>>>(
            relation, std::forward<Args>(args)...);
    };

    template <template <typename...> class InteractionType, typename... ControlParameters,
              template <typename...> class RelationType, typename... RelationParameters, typename... Args>
    auto &addInteractionDynamicsWithInitialization(
        RelationType<RelationParameters...> &relation, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            InteractionDynamicsCK<
                ExecutionPolicy, InteractionType<RelationType<WithInitialization, ControlParameters..., RelationParameters...>>>>(
            relation, std::forward<Args>(args)...);
    };

    template <template <typename...> class InteractionType, typename... ControlParameters,
              template <typename...> class RelationType, typename... RelationParameters, typename... Args>
    auto &addRK2Sequence(
        RelationType<RelationParameters...> &relation, Args &&...args)
    {
        return *particle_dynamics_keeper_.createPtr<
            RungeKuttaSequence<InteractionDynamicsCK<
                ExecutionPolicy,
                InteractionType<RelationType<OneLevel, RungeKutta1stStage, ControlParameters..., RelationParameters...>>,
                InteractionType<RelationType<OneLevel, RungeKutta2ndStage, ControlParameters..., RelationParameters...>>>>>(
            relation, std::forward<Args>(args)...);
    };

    template <template <typename...> class RecorderType, typename... Args>
    auto &addBodyStateRecorder(Args &&...args)
    {
        return *state_recorders_keeper_.createPtr<RecorderType<ExecutionPolicy>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class RegressionType, typename... Parameters, typename... Args>
    auto &addObserveRegression(Args &&...args)
    {
        return *other_io_keeper_.createPtr<
            RegressionType<ObservedQuantityRecording<ExecutionPolicy, Parameters...>>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class RegressionType, template <typename...> class LocalReduceMethodType,
              typename... Parameters, class DynamicsIdentifier, typename... Args>
    auto &addReduceRegression(DynamicsIdentifier &dynamics_identifier, Args &&...args)
    {
        return *other_io_keeper_.createPtr<
            RegressionType<ReducedQuantityRecording<
                ExecutionPolicy, LocalReduceMethodType<Parameters..., DynamicsIdentifier>>>>(
            dynamics_identifier, std::forward<Args>(args)...);
    };

    template <template <typename...> class RegressionType, class LocalReduceMethodType, typename... Args>
    auto &addReduceRegression(Args &&...args)
    {
        return *other_io_keeper_.createPtr<
            RegressionType<ReducedQuantityRecording<ExecutionPolicy, LocalReduceMethodType>>>(std::forward<Args>(args)...);
    };

    template <template <typename...> class IOType, typename... Parameters, typename... Args>
    auto &addIODynamics(Args &&...args)
    {
        return *other_io_keeper_.createPtr<
            IOType<ExecutionPolicy, Parameters...>>(std::forward<Args>(args)...);
    };
};
} // namespace SPH
#endif // PARTICLE_METHOD_CONTAINER_H