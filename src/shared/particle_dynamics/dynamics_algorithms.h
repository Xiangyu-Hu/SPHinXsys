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
 * @file 	dynamics_algorithms.h
 * @brief 	This is the classes for algorithms particle dynamics .
 * @detail	Generally, there are two types of particle dynamics algorithms.
 *			One leads to the change of particle states, the other not.
 *			There are 5 classes the first type. They are:
 * 			SimpleDynamics is without particle interaction. Particles just update their states;
 *			InteractionDynamics is with particle interaction with its neighbors;
 *			InteractionSplit is InteractionDynamics but using spliting algorithm;
 *			InteractionWithUpdate is with particle interaction with its neighbors and then update their states;
 *			Dynamics1Level is the most complex dynamics, has successive three steps: initialization, interaction and update.
 *			In order to avoid misusing of the above algorithms, type traits are used to make sure that the matching between
 *			the algorithm and local dynamics. For example, the LocalDynamics which matches InteractionDynamics must have
 *			the function interaction() but should not have the function update() or initialize().
 *			The existence of the latter suggests that more complex algorithms,
 *			such as InteractionWithUpdate or Dynamics1Level should be used.
 *			There are 2 classes for the second type.
 *			ReduceDynamics carries out a reduce operation through the particles.
 *			Average further computes average of a ReduceDynamics for summation.
 *			Each particle dynamics is templated with a LocalDynamics and a DynamicsRange.
 *			The local dynamics defines the behavior of a single particle or with its neighbors,
 *			and is recognized by particle dynamics with the signature functions, like update, initialization and interaction.
 *			DynamicsRange define and range of particles for the dynamics.
 *			The default range is the entire body. Other ranges are BodyPartByParticle and BodyPartByCell.
 * @author	Chi Zhang, Fabien Pean and Xiangyu Hu
 */

#ifndef DYNAMICS_ALGORITHMS_H
#define DYNAMICS_ALGORITHMS_H

#include "base_dynamics.h"
#include "base_local_dynamics.h"
#include "base_particle_dynamics.h"
#include "cell_linked_list.hpp"
#include "particle_iterators.h"

#include <type_traits>

namespace SPH
{
template <class T, class = void>
struct has_initialize : std::false_type
{
};

template <class T>
struct has_initialize<T, std::void_t<decltype(&T::initialize)>> : std::true_type
{
};

template <class T, class = void>
struct has_interaction : std::false_type
{
};

template <class T>
struct has_interaction<T, std::void_t<decltype(&T::interaction)>> : std::true_type
{
};

template <class T, class = void>
struct has_update : std::false_type
{
};

template <class T>
struct has_update<T, std::void_t<decltype(&T::update)>> : std::true_type
{
};

using namespace execution;

/**
 * @class SimpleDynamics
 * @brief Simple particle dynamics without considering particle interaction
 */
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class SimpleDynamics : public LocalDynamicsType, public BaseDynamics<void>
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    SimpleDynamics(DynamicsIdentifier &identifier, Args &&...args)
        : LocalDynamicsType(identifier, std::forward<Args>(args)...),
          BaseDynamics<void>()
    {
        static_assert(!has_initialize<LocalDynamicsType>::value &&
                          !has_interaction<LocalDynamicsType>::value,
                      "LocalDynamicsType does not fulfill SimpleDynamics requirements");
    };
    virtual ~SimpleDynamics() {};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_->getSPHBody());
        this->setupDynamics(dt);
        particle_for(ExecutionPolicy(),
                     this->identifier_->LoopRange(),
                     [&](size_t i)
                     { this->update(i, dt); });
    };
};

/**
 * @class ReduceDynamics
 * @brief Template class for particle-wise reduce operation, summation, max or min.
 */
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class ReduceDynamics : public LocalDynamicsType,
                       public BaseDynamics<typename LocalDynamicsType::ReturnType>
{
  public:
    using ReturnType = typename LocalDynamicsType::ReturnType;
    template <class DynamicsIdentifier, typename... Args>
    ReduceDynamics(DynamicsIdentifier &identifier, Args &&...args)
        : LocalDynamicsType(identifier, std::forward<Args>(args)...),
          BaseDynamics<ReturnType>(){};
    virtual ~ReduceDynamics() {};
    std::string QuantityName() { return this->quantity_name_; };

    virtual ReturnType exec(Real dt = 0.0) override
    {
        this->setupDynamics(dt);
        ReturnType temp = particle_reduce(ExecutionPolicy(),
                                          this->identifier_->LoopRange(), this->Reference(), this->getOperation(),
                                          [&](size_t i) -> ReturnType
                                          { return this->reduce(i, dt); });
        return this->outputResult(temp);
    };
};

/**
 * @class BaseInteractionDynamics
 * @brief This is the base class for particle interaction with other particles
 */
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class BaseInteractionDynamics : public LocalDynamicsType, public BaseDynamics<void>
{
  public:
    template <typename... Args>
    explicit BaseInteractionDynamics(Args &&...args)
        : LocalDynamicsType(std::forward<Args>(args)...),
          BaseDynamics<void>(){};
    virtual ~BaseInteractionDynamics() {};

    /** pre process such as update ghost state */
    StdVec<BaseDynamics<void> *> pre_processes_;
    /** post process such as impose constraint */
    StdVec<BaseDynamics<void> *> post_processes_;
    /** run the main interaction step between particles. */
    virtual void runMainStep(Real dt) = 0;

    /** run the interactions between particles. */
    virtual void runInteraction(Real dt)
    {
        for (size_t k = 0; k < this->pre_processes_.size(); ++k)
            this->pre_processes_[k]->exec(dt);

        runMainStep(dt);

        for (size_t k = 0; k < this->post_processes_.size(); ++k)
            this->post_processes_[k]->exec(dt);
    };

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_->getSPHBody());
        this->setupDynamics(dt);
        runInteraction(dt);
    };
};

/**
 * @class InteractionSplit
 * @brief This is for the splitting algorithm
 */
template <class LocalDynamicsType, class CellLinkedListType, class ExecutionPolicy = ParallelPolicy>
class BaseInteractionSplit : public BaseInteractionDynamics<LocalDynamicsType, ParallelPolicy>
{
  protected:
    RealBody &real_body_;
    CellLinkedListType &cell_linked_list_;

  public:
    template <typename... Args>
    explicit BaseInteractionSplit(Args &&...args)
        : BaseInteractionDynamics<LocalDynamicsType, ParallelPolicy>(std::forward<Args>(args)...),
          real_body_(DynamicCast<RealBody>(this, this->getSPHBody())),
          cell_linked_list_(DynamicCast<CellLinkedListType>(this, real_body_.getCellLinkedList()))
    {
        static_assert(!has_initialize<LocalDynamicsType>::value &&
                          !has_update<LocalDynamicsType>::value,
                      "LocalDynamicsType does not fulfill InteractionSplit requirements");
    };

    /** run the main interaction step between particles. */
    void runMainStep(Real dt) override
    {
        cell_linked_list_.particle_for_split(ExecutionPolicy(), [&](size_t i)
                                             { this->interaction(i, dt * 0.5); });
    }
};

template <class LocalDynamicsType>
using InteractionSplit = BaseInteractionSplit<LocalDynamicsType, CellLinkedList>;

template <class LocalDynamicsType>
using InteractionAdaptiveSplit = BaseInteractionSplit<LocalDynamicsType, MultilevelCellLinkedList>;

/**
 * @class InteractionDynamics
 * @brief This is the class with a single step of particle interaction with other particles
 */
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class InteractionDynamics : public BaseInteractionDynamics<LocalDynamicsType, ExecutionPolicy>
{
  public:
    template <typename... Args>
    InteractionDynamics(Args &&...args)
        : InteractionDynamics(true, std::forward<Args>(args)...)
    {
        static_assert(!has_initialize<LocalDynamicsType>::value &&
                          !has_update<LocalDynamicsType>::value,
                      "LocalDynamicsType does not fulfill InteractionDynamics requirements");
    };
    virtual ~InteractionDynamics() {};

    /** run the main interaction step between particles. */
    virtual void runMainStep(Real dt) override
    {
        particle_for(ExecutionPolicy(),
                     this->identifier_->LoopRange(),
                     [&](size_t i)
                     { this->interaction(i, dt); });
    }

  protected:
    template <typename... Args>
    InteractionDynamics(bool mostDerived, Args &&...args)
        : BaseInteractionDynamics<LocalDynamicsType, ExecutionPolicy>(std::forward<Args>(args)...){};
};

/**
 * @class InteractionWithUpdate
 * @brief This class includes an interaction and a update steps
 */
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class InteractionWithUpdate : public InteractionDynamics<LocalDynamicsType, ExecutionPolicy>
{
  public:
    template <typename... Args>
    InteractionWithUpdate(Args &&...args)
        : InteractionDynamics<LocalDynamicsType, ExecutionPolicy>(false, std::forward<Args>(args)...)
    {
        static_assert(!has_initialize<LocalDynamicsType>::value,
                      "LocalDynamicsType does not fulfill InteractionWithUpdate requirements");
    }
    virtual ~InteractionWithUpdate() {};

    virtual void exec(Real dt = 0.0) override
    {
        InteractionDynamics<LocalDynamicsType, ExecutionPolicy>::exec(dt);
        particle_for(ExecutionPolicy(),
                     this->identifier_->LoopRange(),
                     [&](size_t i)
                     { this->update(i, dt); });
    };
};

/**
 * @class InteractionWithInitialization
 * @brief This class includes an interaction and an initialization steps
 */
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class InteractionWithInitialization : public InteractionDynamics<LocalDynamicsType, ExecutionPolicy>
{
  public:
    template <typename... Args>
    InteractionWithInitialization(Args &&...args)
        : InteractionDynamics<LocalDynamicsType, ExecutionPolicy>(false, std::forward<Args>(args)...)
    {
        static_assert(!has_update<LocalDynamicsType>::value,
                      "LocalDynamicsType does not fulfill InteractionWithInitialization requirements");
    }
    virtual ~InteractionWithInitialization() {};

    virtual void exec(Real dt = 0.0) override
    {
        particle_for(ExecutionPolicy(),
                     this->identifier_->LoopRange(),
                     [&](size_t i)
                     { this->initialization(i, dt); });
        InteractionDynamics<LocalDynamicsType, ExecutionPolicy>::exec(dt);
    };
};

/**
 * @class Dynamics1Level
 * @brief This class includes three steps, including initialization, interaction and update.
 * It is the most complex particle dynamics type,
 * and is typically for computing the main fluid and solid dynamics.
 */
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class Dynamics1Level : public InteractionDynamics<LocalDynamicsType, ExecutionPolicy>
{
  public:
    template <typename... Args>
    Dynamics1Level(Args &&...args)
        : InteractionDynamics<LocalDynamicsType, ExecutionPolicy>(
              false, std::forward<Args>(args)...) {}
    virtual ~Dynamics1Level() {};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_->getSPHBody());
        this->setupDynamics(dt);

        particle_for(ExecutionPolicy(),
                     this->identifier_->LoopRange(),
                     [&](size_t i)
                     { this->initialization(i, dt); });

        InteractionDynamics<LocalDynamicsType, ExecutionPolicy>::runInteraction(dt);

        particle_for(ExecutionPolicy(),
                     this->identifier_->LoopRange(),
                     [&](size_t i)
                     { this->update(i, dt); });
    };
};
} // namespace SPH
#endif // DYNAMICS_ALGORITHMS_H
