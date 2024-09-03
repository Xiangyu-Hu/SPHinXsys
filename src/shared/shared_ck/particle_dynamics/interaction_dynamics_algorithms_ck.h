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
class InteractionDynamicsCK<Base>
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
};

template <>
class InteractionDynamicsCK<WithUpdate> : public InteractionDynamicsCK<Base>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<Base>(){};

    virtual void runInteraction(Real dt) override
    {
        InteractionDynamicsCK<Base>::runInteraction(dt);
        runUpdateStep(dt);
    };

  protected:
    virtual void runUpdateStep(Real dt) = 0;
};

template <>
class InteractionDynamicsCK<WithInitialization> : public InteractionDynamicsCK<Base>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<Base>(){};

    virtual void runInteraction(Real dt) override
    {
        runInitialization(dt);
        InteractionDynamicsCK<Base>::runInteraction(dt);
    };

  protected:
    virtual void runInitialization(Real dt) = 0;
};

template <>
class InteractionDynamicsCK<OneLevel> : public InteractionDynamicsCK<Base>
{
  public:
    InteractionDynamicsCK() : InteractionDynamicsCK<Base>(){};

    virtual void runInteraction(Real dt) override
    {
        runInitialization(dt);
        InteractionDynamicsCK<Base>::runInteraction(dt);
        runUpdateStep(dt);
    };

  protected:
    virtual void runUpdateStep(Real dt) = 0;
    virtual void runInitialization(Real dt) = 0;
};

template <class ExecutionPolicy, class InteractionType,
          template <typename...> class InteractionName, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, InteractionType,
                            InteractionName<Inner<Parameters...>>>
    : public InteractionName<Inner<Parameters...>>,
      public InteractionDynamicsCK<InteractionType>
{

    using LocalDynamicsType = InteractionName<Inner<Parameters...>>;
    using ComputingKernel = typename LocalDynamicsType::ComputingKernel;
    using KernelImplementation = Implementation<LocalDynamicsType, ExecutionPolicy>;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionName<Inner<Parameters...>>(std::forward<Args>(args)...),
          InteractionDynamicsCK<InteractionType>(),
          ex_policy_(ExecutionPolicy{}), kernel_implementation_(*this){};
    virtual ~InteractionDynamicsCK(){};

  protected:
    ExecutionPolicy ex_policy_;
    KernelImplementation kernel_implementation_;

    ComputingKernel *getFirstComputingKernel()
    {
        return kernel_implementation_->getComputingKernel();
    };

    virtual void runMainStep(Real dt) override
    {
        ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();
        particle_for(ex_policy_,
                     this->identifier_.LoopRange(),
                     [=](size_t i)
                     { computing_kernel->interaction(i, dt); });
    };
};

template <class ExecutionPolicy, class InteractionType,
          template <typename...> class InteractionName, typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy, InteractionType,
                            InteractionName<Contact<Parameters...>>>
    : public InteractionName<Contact<Parameters...>>,
      public InteractionDynamicsCK<InteractionType>
{

    using LocalDynamicsType = InteractionName<Contact<Parameters...>>;
    using ComputingKernel = typename LocalDynamicsType::ComputingKernel;
    using KernelImplementation = Implementation<LocalDynamicsType, ExecutionPolicy>;

  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionName<Contact<Parameters...>>(std::forward<Args>(args)...),
          InteractionDynamicsCK<InteractionType>(),
          ex_policy_(ExecutionPolicy{})
    {
        for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
        {
            contact_kernel_implementation_.push_back(
                contact_kernel_implementation_ptrs_.template createPtr<KernelImplementation>(*this));
        }
    };
    virtual ~InteractionDynamicsCK(){};

  protected:
    ExecutionPolicy ex_policy_;
    StdVec<KernelImplementation *> contact_kernel_implementation_;
    ComputingKernel *getFirstComputingKernel()
    {
        return contact_kernel_implementation_[0]->getComputingKernel(0);
    };

    virtual void runMainStep(Real dt) override
    {
        for (size_t k = 0; k != this->contact_bodies_.size(); ++k)
        {
            ComputingKernel *computing_kernel = contact_kernel_implementation_[k]->getComputingKernel(k);

            particle_for(ex_policy_,
                         this->identifier_.LoopRange(),
                         [=](size_t i)
                         { computing_kernel->interaction(i, dt); });
        }
    };
};

template <class ExecutionPolicy,
          template <typename...> class InteractionName,
          template <typename...> class RelationType,
          typename... Parameters>
class InteractionDynamicsCK<ExecutionPolicy,
                            InteractionName<RelationType<Parameters...>>>
    : public InteractionDynamicsCK<ExecutionPolicy, Base,
                                   InteractionName<RelationType<Parameters...>>>,
      public BaseDynamics<void>
{
  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionDynamicsCK<ExecutionPolicy, Base,
                                InteractionName<RelationType<Parameters...>>>(std::forward<Args>(args)...),
          BaseDynamics<void>() {}
    virtual ~InteractionDynamicsCK(){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_.getSPHBody());
        this->setupDynamics(dt);
        runInteraction(dt);
    };
};

template <class ExecutionPolicy,
          template <typename...> class InteractionName,
          template <typename...> class RelationType,
          typename... OtherParameters>
class InteractionDynamicsCK<ExecutionPolicy,
                            InteractionName<RelationType<WithUpdate, OtherParameters...>>>
    : public InteractionDynamicsCK<ExecutionPolicy, WithUpdate,
                                   InteractionName<RelationType<WithUpdate, OtherParameters...>>>,
      public BaseDynamics<void>
{
  public:
    template <typename... Args>
    InteractionDynamicsCK(Args &&...args)
        : InteractionDynamicsCK<
              ExecutionPolicy, WithUpdate,
              InteractionName<RelationType<WithUpdate, OtherParameters...>>>(std::forward<Args>(args)...),
          BaseDynamics<void>() {}
    virtual ~InteractionDynamicsCK(){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_.getSPHBody());
        this->setupDynamics(dt);
        runInteraction(dt);
    };

  protected:
    virtual void runUpdateStep(Real dt) override
    {
        ComputingKernel *computing_kernel = this->getFirstComputingKernel();
        particle_for(this->ex_policy_,
                     this->identifier_.LoopRange(),
                     [=](size_t i)
                     { computing_kernel->update(i, dt); });
    };
};
} // namespace SPH
#endif // INTERACTION_DYNAMICS_ALGORITHMS_CK_H
