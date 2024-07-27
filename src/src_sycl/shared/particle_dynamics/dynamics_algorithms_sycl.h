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
 * @file 	dynamics_algorithms_sycl.h
 * @brief TBD.
 * @author Alberto Guarnieri and Xiangyu Hu
 */

#ifndef DYNAMICS_ALGORITHMS_SYCL_H
#define DYNAMICS_ALGORITHMS_SYCL_H

#include "dynamics_algorithms.h"
#include "particle_iterators_sycl.h"

namespace SPH
{
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class SimpleDynamicsSYCL : public LocalDynamicsType, public BaseDynamics<void>
{
    using ComputingKernelType = LocalDynamicsType::ComputingKernelType;

  public:
    template <class DynamicsIdentifier, typename... Args>
    SimpleDynamicsSYCL(DynamicsIdentifier &identifier, Args &&...args)
        : LocalDynamicsType(ExecutionPolicy{}, identifier, std::forward<Args>(args)...),
          BaseDynamics<void>(identifier.getSPHBody()),
          computing_kernel_implementation_(this->computing_kernel_)
    {
        static_assert(!has_initialize<ComputingKernelType>::value &&
                          !has_interaction<ComputingKernelType>::value,
                      "LocalDynamicsType does not fulfill SimpleDynamics requirements");
    };
    virtual ~SimpleDynamicsSYCL(){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated();
        this->setupDynamics(dt);
        particle_for(computing_kernel_implementation_,
                     this->LoopRange(),
                     [&](size_t i, auto &&computing_kernel)
                     { computing_kernel.update(i, dt); });
    };

  protected:
    Implementation<ComputingKernelType, ExecutionPolicy> computing_kernel_implementation_;
};
} // namespace SPH
#endif // DYNAMICS_ALGORITHMS_SYCL_H
