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
 * @file 	simple_algorithms_ck.h
 * @brief 	This is the classes for algorithms particle dynamics .
 * @detail	TBD
 * @author	Xiangyu Hu
 */

#ifndef SIMPLE_ALGORITHMS_CK_H
#define SIMPLE_ALGORITHMS_CK_H

#include "base_local_dynamics.h"
#include "base_particle_dynamics.h"
#include "particle_iterators_ck.h"

namespace SPH
{
template <class ExecutionPolicy, class UpdateType>
class StateDynamics : public UpdateType, public BaseDynamics<void>
{
    using UpdateKernel = typename UpdateType::UpdateKernel;
    using KernelImplementation =
        Implementation<ExecutionPolicy, UpdateType, UpdateKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <class DynamicsIdentifier, typename... Args>
    StateDynamics(DynamicsIdentifier &identifier, Args &&...args)
        : UpdateType(identifier, std::forward<Args>(args)...),
          BaseDynamics<void>(), kernel_implementation_(*this){};
    virtual ~StateDynamics(){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_.getSPHBody());
        this->setupDynamics(dt);
        auto &loop_range = kernel_implementation_.getLoopRange();
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        particle_for(loop_range, [=](size_t i)
                     { update_kernel->update(i, dt); });
    };
};

template <class ExecutionPolicy, class ReduceType>
class ReduceDynamicsCK : public ReduceType,
                         public BaseDynamics<typename ReduceType::ReturnType>
{
    using ReduceKernel = typename ReduceType::ReduceKernel;
    using ReturnType = typename ReduceType::ReturnType;
    using KernelImplementation =
        Implementation<ExecutionPolicy, ReduceType, ReduceKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <class DynamicsIdentifier, typename... Args>
    ReduceDynamicsCK(DynamicsIdentifier &identifier, Args &&...args)
        : ReduceType(identifier, std::forward<Args>(args)...),
          BaseDynamics<ReturnType>(), kernel_implementation_(*this){};
    virtual ~ReduceDynamicsCK(){};

    std::string QuantityName() { return this->quantity_name_; };
    std::string DynamicsIdentifierName() { return this->identifier_.getName(); };

    virtual ReturnType exec(Real dt = 0.0) override
    {
        this->setupDynamics(dt);
        auto &loop_range = kernel_implementation_.getLoopRange();
        ReduceKernel *reduce_kernel = kernel_implementation_.getComputingKernel();
        ReturnType temp = particle_reduce(
            loop_range, this->Reference(), this->getOperation(),
            [=](size_t i) -> ReturnType
            { return reduce_kernel->reduce(i, dt); });
        return this->outputResult(temp);
    };
};
} // namespace SPH
#endif // SIMPLE_ALGORITHMS_CK_H
