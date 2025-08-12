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
 * @file kernel_gradient_integral.h
 * @brief The error residual for computing the gradient of unit.
 * @author	Xiangyu Hu
 */

#ifndef ZERO_GRADIENT_RESIDUAL_H
#define ZERO_GRADIENT_RESIDUAL_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class BaseInteractionType>
class KernelGradientIntegralBase : public BaseInteractionType
{
  public:
    template <class DynamicsIdentifier>
    explicit KernelGradientIntegralBase(DynamicsIdentifier &identifier);
    virtual ~KernelGradientIntegralBase() {}

  protected:
    DiscreteVariable<Vecd> *dv_kernel_gradient_integral_; ///< "KernelGradientIntegral"
};

template <typename...>
class KernelGradientIntegral;

template <class KernelCorrectionType, typename... Parameters>
class KernelGradientIntegral<Inner<KernelCorrectionType, Parameters...>>
    : public KernelGradientIntegralBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = KernelGradientIntegralBase<Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit KernelGradientIntegral(Inner<Parameters...> &inner_relation);
    virtual ~KernelGradientIntegral() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Vecd *kernel_gradient_integral_;
        Real *Vol_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
};

template <class KernelCorrectionType, typename... Parameters>
class KernelGradientIntegral<Contact<Boundary, KernelCorrectionType, Parameters...>>
    : public KernelGradientIntegralBase<Interaction<Contact<Parameters...>>>
{
    using BaseInteraction = KernelGradientIntegralBase<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit KernelGradientIntegral(Contact<Parameters...> &contact_relation);
    virtual ~KernelGradientIntegral() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Vecd *kernel_gradient_integral_;
        Real *contact_Vol_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
};
} // namespace SPH
#endif // ZERO_GRADIENT_RESIDUAL_H
