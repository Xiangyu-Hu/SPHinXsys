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
 * @file 	interaction_algorithms_ck.h
 * @brief 	TBD
 * @author	Xiangyu Hu
 */

#ifndef INTERACTION_ALGORITHMS_CK_H
#define INTERACTION_ALGORITHMS_CK_H

#include "base_particle_dynamics.h"
#include "interaction_ck.hpp"
#include "particle_iterators_ck.h"
#include "simple_algorithms_ck.h"

namespace SPH
{
template <typename...>
class InteractionDynamicsCK;

template <>
class InteractionDynamicsCK<Base>
{
  public:
    InteractionDynamicsCK() {};

  protected:
    /** pre process such as update ghost state */
    StdVec<BaseDynamics<void> *> pre_processes_;
    /** post process such as impose constraint */
    StdVec<BaseDynamics<void> *> post_processes_;
    /** run the interaction step between particles. */
    virtual void runInteractionStep(Real dt) = 0;
    /** run all interactions step. */
    virtual void runAllSteps(Real dt);
};

template <>
class InteractionDynamicsCK<WithUpdate> : public InteractionDynamicsCK<Base>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<Base>() {};
    virtual void runAllSteps(Real dt) override;

  protected:
    virtual void runUpdateStep(Real dt) = 0;
};

template <>
class InteractionDynamicsCK<WithInitialization> : public InteractionDynamicsCK<Base>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<Base>() {};
    virtual void runAllSteps(Real dt) override;

  protected:
    virtual void runInitializationStep(Real dt) = 0;
};

template <>
class InteractionDynamicsCK<OneLevel> : public InteractionDynamicsCK<Base>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<Base>() {};
    virtual void runAllSteps(Real dt) override;

  protected:
    virtual void runInitializationStep(Real dt) = 0;
    virtual void runUpdateStep(Real dt) = 0;
};

template <class ExecutionPolicy, typename AlgorithmType, template <typename...> class InteractionType>
class InteractionDynamicsCK<ExecutionPolicy, InteractionType<AlgorithmType>>
    : public InteractionDynamicsCK<AlgorithmType>,
      public BaseDynamics<void>
{
    UniquePtrsKeeper<BaseDynamics<void>> supplementary_dynamics_keeper_;

  public:
    InteractionDynamicsCK() {};

    template <typename... ControlParameters, typename... RelationParameters, typename... Args>
    auto &addPostContactInteraction(Contact<RelationParameters...> &contact_relation, Args &&...args);

    template <typename... ControlParameters, typename... RelationParameters, typename... Args>
    auto &addPreContactInteraction(Contact<RelationParameters...> &contact_relation, Args &&...args);

    auto &addPostContactInteraction(BaseDynamics<void> &contact_interaction);
    auto &addPreContactInteraction(BaseDynamics<void> &contact_interaction);

    template <class UpdateType, typename... Args>
    auto &addPostStateDynamics(Args &&...args);
    auto &addPostStateDynamics(BaseDynamics<void> &state_dynamics);

    template <class UpdateType, typename... Args>
    auto &addPreStateDynamics(Args &&...args);
    auto &addPreStateDynamics(BaseDynamics<void> &state_dynamics);
};

template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, Base, InteractionType<Inner<Parameters...>>>
    : public InteractionType<Inner<Parameters...>>
{
    using LocalDynamicsType = InteractionType<Inner<Parameters...>>;
    using Identifier = typename LocalDynamicsType::Identifier;
    using InteractKernel = typename LocalDynamicsType::InteractKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args);
    virtual ~InteractionDynamicsCK() {};

  protected:
    void runInteraction(Real dt);
};

template <class ExecutionPolicy, template <typename...> class InteractionType, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, Base, InteractionType<Contact<Parameters...>>>
    : public InteractionType<Contact<Parameters...>>
{
    using LocalDynamicsType = InteractionType<Contact<Parameters...>>;
    using Identifier = typename LocalDynamicsType::Identifier;
    using InteractKernel = typename LocalDynamicsType::InteractKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;
    UniquePtrsKeeper<KernelImplementation> contact_kernel_implementation_ptrs_;
    StdVec<KernelImplementation *> contact_kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args);
    virtual ~InteractionDynamicsCK() {};

  protected:
    void runInteraction(Real dt);
};

template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, InteractionType<RelationType<Parameters...>>>
    : public InteractionDynamicsCK<
          ExecutionPolicy, Base, InteractionType<RelationType<Parameters...>>>,
      public InteractionDynamicsCK<ExecutionPolicy, InteractionType<Base>>
{
  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args);
    virtual ~InteractionDynamicsCK() {};
    virtual void exec(Real dt = 0.0) override;
    virtual void runInteractionStep(Real dt = 0.0) override;
};

template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
class InteractionDynamicsCK<
    ExecutionPolicy, InteractionType<RelationType<WithUpdate, OtherParameters...>>>
    : public InteractionDynamicsCK<
          ExecutionPolicy, Base, InteractionType<RelationType<WithUpdate, OtherParameters...>>>,
      public InteractionDynamicsCK<ExecutionPolicy, InteractionType<WithUpdate>>
{
    using LocalDynamicsType = InteractionType<RelationType<WithUpdate, OtherParameters...>>;
    using Identifier = typename LocalDynamicsType::Identifier;
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using BaseInteractKernel = typename LocalDynamicsType::BaseInteractKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args);
    virtual ~InteractionDynamicsCK() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    virtual void runInteractionStep(Real dt = 0.0) override;
    virtual void runUpdateStep(Real dt) override;
};

template <class ExecutionPolicy, template <typename...> class InteractionType,
          template <typename...> class RelationType, typename... OtherParameters>
class InteractionDynamicsCK<
    ExecutionPolicy, InteractionType<RelationType<OneLevel, OtherParameters...>>>
    : public InteractionDynamicsCK<
          ExecutionPolicy, Base, InteractionType<RelationType<OneLevel, OtherParameters...>>>,
      public InteractionDynamicsCK<ExecutionPolicy, InteractionType<OneLevel>>
{
    using LocalDynamicsType = InteractionType<RelationType<OneLevel, OtherParameters...>>;
    using Identifier = typename LocalDynamicsType::Identifier;
    using InitializeKernel = typename LocalDynamicsType::InitializeKernel;
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using BaseInteractKernel = typename LocalDynamicsType::BaseInteractKernel;
    using InitializeKernelImplementation =
        Implementation<ExecutionPolicy, LocalDynamicsType, InitializeKernel>;
    using UpdateKernelImplementation =
        Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;

    InitializeKernelImplementation initialize_kernel_implementation_;
    UpdateKernelImplementation update_kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args);
    virtual ~InteractionDynamicsCK() {};
    virtual void exec(Real dt = 0.0) override;

  protected:
    virtual void runInitializationStep(Real dt) override;
    virtual void runInteractionStep(Real dt = 0.0) override;
    virtual void runUpdateStep(Real dt) override;
};

template <class ExecutionPolicy, template <typename...> class InteractionType>
class InteractionDynamicsCK<ExecutionPolicy, InteractionType<>>
{
  public:
    InteractionDynamicsCK() {};
    void runInteractionStep(Real dt = 0.0) {};
};

template <class ExecutionPolicy, template <typename...> class InteractionType,
          class FirstInteraction, class... Others>
class InteractionDynamicsCK<ExecutionPolicy, InteractionType<FirstInteraction, Others...>>
    : public InteractionDynamicsCK<ExecutionPolicy, InteractionType<FirstInteraction>>
{
  protected:
    InteractionDynamicsCK<ExecutionPolicy, InteractionType<Others...>> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit InteractionDynamicsCK(
        FirstParameterSet &&first_parameter_set, OtherParameterSets &&...other_parameter_sets);
    virtual void runInteractionStep(Real dt = 0.0) override;
};
} // namespace SPH
#endif // INTERACTION_ALGORITHMS_CK_H
