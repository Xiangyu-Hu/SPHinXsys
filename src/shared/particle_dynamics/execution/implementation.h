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
 * @file 	implementation.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef IMPLEMENTATION_H
#define IMPLEMENTATION_H

#include "base_implementation.h"
#include "base_data_type.h"
#include "execution_policy.h"
#include "loop_range.h"
#include "ownership.h"

namespace SPH
{
namespace execution
{
template <class ExecutionPolicy, class LocalDynamicsType>
class Implementation<ExecutionPolicy, LocalDynamicsType> : public Implementation<Base>
{
    using DynamicsIdentifier = typename LocalDynamicsType::DynamicsIdentifierType;
    using LoopRangeType = LoopRangeCK<ExecutionPolicy, DynamicsIdentifier>;

  public:
    explicit Implementation(LocalDynamicsType &local_dynamics)
        : Implementation<Base>(), local_dynamics_(local_dynamics),
          loop_range_(nullptr) {};
    ~Implementation()
    {
        delete loop_range_;
    };

    LoopRangeType &getLoopRange()
    {
        if (loop_range_ == nullptr)
        {
            loop_range_ = new LoopRangeType(local_dynamics_.getDynamicsIdentifier());
        }
        return *loop_range_;
    };

  protected:
    LocalDynamicsType &local_dynamics_;
    LoopRangeType *loop_range_;
};

template <class ExecutionPolicy, class LocalDynamicsType, class ComputingKernelType>
class Implementation<ExecutionPolicy, LocalDynamicsType, ComputingKernelType>
    : public Implementation<ExecutionPolicy, LocalDynamicsType>
{
  public:
    explicit Implementation(LocalDynamicsType &local_dynamics)
        : Implementation<ExecutionPolicy, LocalDynamicsType>(local_dynamics),
          computing_kernel_(nullptr) {}
    ~Implementation()
    {
        delete computing_kernel_;
    }

    template <typename... Args>
    ComputingKernelType *getComputingKernel(Args &&...args)
    {
        if (computing_kernel_ == nullptr)
        {
            this->local_dynamics_.registerComputingKernel(this, std::forward<Args>(args)...);
            computing_kernel_ = new ComputingKernelType(
                ExecutionPolicy{}, this->local_dynamics_, std::forward<Args>(args)...);
            this->setUpdated();
        }

        if (!this->isUpdated())
        {
            overwriteComputingKernel(std::forward<Args>(args)...);
        }

        return computing_kernel_;
    }

    template <typename... Args>
    void overwriteComputingKernel(Args &&...args)
    {
        *computing_kernel_ = ComputingKernelType(
            ExecutionPolicy{}, this->local_dynamics_, std::forward<Args>(args)...);
        this->setUpdated();
    }

  protected:
    ComputingKernelType *computing_kernel_;
};
} // namespace execution
} // namespace SPH
#endif // IMPLEMENTATION_H
