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
 * @file kernel_integral.h
 * @brief The integral of the kernel function throughout the support.
 * @author	Xiangyu Hu
 */

#ifndef KERNEL_INTEGRAL_H
#define KERNEL_INTEGRAL_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class BaseInteractionType>
class KernelIntegralBase : public BaseInteractionType
{
  public:
    template <class DynamicsIdentifier>
    explicit KernelIntegralBase(DynamicsIdentifier &identifier);
    virtual ~KernelIntegralBase() {}

  protected:
    DiscreteVariable<Real> *dv_kernel_integral_; ///< "KernelIntegral"
};

template <typename...>
class KernelIntegral;

template <typename... Parameters>
class KernelIntegral<Inner<Parameters...>>
    : public KernelIntegralBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = KernelIntegralBase<Interaction<Inner<Parameters...>>>;

  public:
    explicit KernelIntegral(Inner<Parameters...> &inner_relation);
    virtual ~KernelIntegral() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real W0_;
        Real *kernel_integral_;
        Real *Vol_;
    };
};

template <typename... Parameters>
class KernelIntegral<Contact<Parameters...>>
    : public KernelIntegralBase<Interaction<Contact<Parameters...>>>
{
    using BaseInteraction = KernelIntegralBase<Interaction<Contact<Parameters...>>>;

  public:
    explicit KernelIntegral(Contact<Parameters...> &contact_relation);
    virtual ~KernelIntegral() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *kernel_integral_;
        Real *contact_Vol_;
    };
};
} // namespace SPH
#endif // KERNEL_INTEGRAL_H
