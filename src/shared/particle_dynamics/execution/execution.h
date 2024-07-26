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
 * @file 	execution.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef EXECUTION_H
#define EXECUTION_H

#include "execution_policy.h"

namespace SPH
{
namespace execution
{
template <typename... T>
class Implementation;

template <class ComputingKernelType, class ExecutionPolicy>
class Implementation<ComputingKernelType, ExecutionPolicy>
{
  public:
    Implementation() = default;

    template <class... Args>
    explicit Implementation(const ExecutionPolicy &execution_policy, Args &&...args)
        : computing_kernel_(std::forward<Args>(args)...) {}

    ComputingKernelType &getBuffer()
    {
        return computing_kernel_;
    }

  private:
    ComputingKernelType computing_kernel_;
};
} // namespace execution
} // namespace SPH
#endif // EXECUTION_H
