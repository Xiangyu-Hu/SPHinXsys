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
 * @file 	interaction_dynamics_algorithms_ck.h
 * @brief 	TBD
 * @author	Xiangyu Hu
 */

#ifndef INTERACTION_DYNAMICS_ALGORITHMS_CK_H
#define INTERACTION_DYNAMICS_ALGORITHMS_CK_H

#include "dynamics_algorithms.h"

namespace SPH
{
template <typename...>
class InteractionDynamicsCK;

template <>
class InteractionDynamicsCK<>
{
  public:
    InteractionDynamicsCK(){};
    void addPreProcess(BaseDynamics<void> *pre_process) { pre_processes_.push_back(pre_process); };
    void addPostProcess(BaseDynamics<void> *post_process) { post_processes_.push_back(post_process); };

  protected:
    /** pre process such as update ghost state */
    StdVec<BaseDynamics<void> *> pre_processes_;
    /** post process such as impose constraint */
    StdVec<BaseDynamics<void> *> post_processes_;
    /** run the interaction step between particles. */
    virtual void runInteractionStep(Real dt) = 0;

    /** run all interactions step. */
    virtual void runAllSteps(Real dt)
    {
        for (size_t k = 0; k < this->pre_processes_.size(); ++k)
            this->pre_processes_[k]->exec(dt);

        runInteractionStep(dt);

        for (size_t k = 0; k < this->post_processes_.size(); ++k)
            this->post_processes_[k]->exec(dt);
    };
};

template <>
class InteractionDynamicsCK<WithUpdate> : public InteractionDynamicsCK<>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<>(){};

    virtual void runAllSteps(Real dt) override
    {
        InteractionDynamicsCK<>::runAllSteps(dt);
        runUpdateStep(dt);
    };

  protected:
    virtual void runUpdateStep(Real dt) = 0;
};

template <>
class InteractionDynamicsCK<WithInitialization> : public InteractionDynamicsCK<>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<>(){};

    virtual void runAllSteps(Real dt) override
    {
        runInitializationStep(dt);
        InteractionDynamicsCK<>::runAllSteps(dt);
    };

  protected:
    virtual void runInitializationStep(Real dt) = 0;
};

template <>
class InteractionDynamicsCK<OneLevel> : public InteractionDynamicsCK<>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<>(){};

    virtual void runAllSteps(Real dt) override
    {
        runInitializationStep(dt);
        InteractionDynamicsCK<>::runAllSteps(dt);
        runUpdateStep(dt);
    };

  protected:
    virtual void runInitializationStep(Real dt) = 0;
    virtual void runUpdateStep(Real dt) = 0;
};

template <class ExecutionPolicy,
          template <typename...> class InteractionType, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, Base,
                            InteractionType<Inner<Parameters...>>>
    : public InteractionType<Inner<Parameters...>>
{
    using LocalDynamicsType = InteractionType<Inner<Parameters...>>;
    using InteractKernel = typename LocalDynamicsType::InteractKernel;
    using KernelImplementation =
        Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionType<Inner<ExtensionTypes..., Parameters...>>(std::forward<Args>(args)...),
          kernel_implementation_(*this){};
    virtual ~InteractionDynamicsCK(){};

  protected:
    void runInteraction(Real dt)
    {
        InteractKernel *interact_kernel = kernel_implementation_.getComputingKernel();
        particle_for(ExecutionPolicy{},
                     this->identifier_.LoopRange(),
                     [=](size_t i)
                     { interact_kernel->interact(i, dt); });
    };
};

template <class ExecutionPolicy,
          template <typename...> class InteractionType, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, Base,
                            InteractionType<Contact<Parameters...>>>
    : public InteractionType<Contact<Parameters...>>
{
    using LocalDynamicsType = InteractionType<Contact<Parameters...>>;
    using InteractKernel = typename LocalDynamicsType::InteractKernel;
    using KernelImplementation =
        Implementation<ExecutionPolicy, LocalDynamicsType, InteractKernel>;
    StdVec<KernelImplementation *> contact_kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionType<Contact<Parameters...>>(std::forward<Args>(args)...),
          InteractionDynamicsCK<FirstParameter>()
    {
        for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
        {
            contact_kernel_implementation_.push_back(
                contact_kernel_implementation_ptrs_
                    .template createPtr<KernelImplementation>(*this));
        }
    };
    virtual ~InteractionDynamicsCK(){};

  protected:
    void runInteraction(Real dt)
    {
        for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
        {
            InteractKernel *interact_kernel =
                contact_kernel_implementation_[k]->getComputingKernel(k);

            particle_for(ExecutionPolicy{},
                         this->identifier_.LoopRange(),
                         [=](size_t i)
                         { interact_kernel->interact(i, dt); });
        }
    };
};

template <class ExecutionPolicy,
          template <typename...> class InteractionType,
          template <typename...> class RelationType,
          typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy,
                            InteractionType<RelationType<Parameters...>>>
    : public InteractionDynamicsCK<ExecutionPolicy, Base,
                                   InteractionType<RelationType<Parameters...>>>,
      public InteractionDynamicsCK<>,
      public BaseDynamics<void>
{
  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionDynamicsCK<
              ExecutionPolicy, Base,
              InteractionType<RelationType<Parameters...>>>(std::forward<Args>(args)...),
          InteractionDynamicsCK<>(), BaseDynamics<void>() {}
    virtual ~InteractionDynamicsCK(){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_.getSPHBody());
        this->setupDynamics(dt);
        InteractionDynamicsCK<>::runAllSteps(dt);
    };

    virtual void runInteractionStep(Real dt = 0.0) override
    {
        this->runInteraction(dt);
    };
};

template <class ExecutionPolicy,
          template <typename...> class InteractionType,
          template <typename...> class RelationType,
          typename... OtherParameters>
class InteractionDynamicsCK<
    ExecutionPolicy,
    InteractionType<RelationType<WithUpdate, OtherParameters...>>>
    : public InteractionDynamicsCK<
          ExecutionPolicy, Base,
          InteractionType<RelationType<WithUpdate, OtherParameters...>>>,
      public InteractionDynamicsCK<WithUpdate>(),
      public BaseDynamics<void>
{
    using LocalDynamicsType = InteractionType<RelationType<WithUpdate, OtherParameters...>>;
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation =
        Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionDynamicsCK<
              ExecutionPolicy, WithUpdate,
              InteractionType<RelationType<WithUpdate, OtherParameters...>>>(
              std::forward<Args>(args)...),
          InteractionDynamicsCK<WithUpdate>(),
          BaseDynamics<void>(),
          kernel_implementation_(*this) {}
    virtual ~InteractionDynamicsCK(){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_.getSPHBody());
        this->setupDynamics(dt);
        InteractionDynamicsCK<WithUpdate>::runAllSteps(dt);
    };

  protected:
    virtual void runInteractionStep(Real dt = 0.0) override
    {
        this->runInteraction(dt);
    };

    virtual void runUpdateStep(Real dt) override
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        particle_for(ExecutionPolicy{},
                     this->identifier_.LoopRange(),
                     [=](size_t i)
                     { update_kernel->update(i, dt); });
    };
};
} // namespace SPH
#endif // INTERACTION_DYNAMICS_ALGORITHMS_CK_H
