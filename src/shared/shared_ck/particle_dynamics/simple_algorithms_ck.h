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
 * @file 	simple_algorithms_ck.h
 * @brief 	This is the classes for algorithms particle dynamics .
 * @detail	TBD
 * @author	Xiangyu Hu
 */

#ifndef SIMPLE_ALGORITHMS_CK_H
#define SIMPLE_ALGORITHMS_CK_H

#include "base_local_dynamics.h"
#include "base_particle_dynamics.h"
#include "io_log.h"
#include "particle_iterators_ck.h"

namespace SPH
{
template <class ExecutionPolicy, class UpdateType>
class StateDynamics : public UpdateType, public BaseDynamics<void>
{
    using Identifier = typename UpdateType::Identifier;
    using UpdateKernel = typename UpdateType::UpdateKernel;
    using FinishDynamics = typename UpdateType::FinishDynamics;
    using KernelImplementation =
        Implementation<ExecutionPolicy, UpdateType, UpdateKernel>;
    KernelImplementation kernel_implementation_;
    FinishDynamics finish_dynamics_;

  public:
    template <typename... Args>
    StateDynamics(Args &&...args)
        : UpdateType(std::forward<Args>(args)...),
          BaseDynamics<void>(), kernel_implementation_(*this), finish_dynamics_(*this){};
    virtual ~StateDynamics() {};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_->getSPHBody());
        this->setupDynamics(dt);
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        particle_for(LoopRangeCK<ExecutionPolicy, Identifier>(*this->identifier_),
                     [=](size_t i)
                     { update_kernel->update(i, dt); });

        finish_dynamics_();

        this->logger_->debug(
            "StateDynamics::exec() for {} at {}",
            type_name<UpdateType>(),
            this->sph_body_->getName());
    };
};

template <class ExecutionPolicy, class ReduceType>
class ReduceDynamicsCK : public ReduceType,
                         public BaseDynamics<typename ReduceType::FinishDynamics::OutputType>
{
    using Identifier = typename ReduceType::Identifier;
    using ReduceKernel = typename ReduceType::ReduceKernel;
    using ReduceReturnType = typename ReduceType::ReturnType;
    using Operation = typename ReduceType::OperationType;
    using FinishDynamics = typename ReduceType::FinishDynamics;
    using KernelImplementation =
        Implementation<ExecutionPolicy, ReduceType, ReduceKernel>;
    KernelImplementation kernel_implementation_;
    ReduceReturnType reduced_value_;
    FinishDynamics finish_dynamics_;

  public:
    using OutputType = typename FinishDynamics::OutputType;

    template <typename... Args>
    ReduceDynamicsCK(Args &&...args)
        : ReduceType(std::forward<Args>(args)...),
          BaseDynamics<OutputType>(), kernel_implementation_(*this),
          reduced_value_(this->reference_), finish_dynamics_(*this){};
    virtual ~ReduceDynamicsCK() {};
    std::string QuantityName() { return this->quantity_name_; };
    ReduceReturnType ReducedValue() { return reduced_value_; };

    virtual OutputType exec(Real dt = 0.0) override
    {
        this->setupDynamics(dt);
        ReduceKernel *reduce_kernel = kernel_implementation_.getComputingKernel();
        reduced_value_ = particle_reduce<Operation>(
            LoopRangeCK<ExecutionPolicy, Identifier>(*this->identifier_),
            this->reference_,
            [=](size_t i)
            { return reduce_kernel->reduce(i, dt); });

        this->logger_->debug(
            "ReduceDynamicsCK::exec() for {} at {}",
            type_name<ReduceType>(),
            this->sph_body_->getName());

        return finish_dynamics_.Result(reduced_value_);
    };
};
} // namespace SPH
#endif // SIMPLE_ALGORITHMS_CK_H
