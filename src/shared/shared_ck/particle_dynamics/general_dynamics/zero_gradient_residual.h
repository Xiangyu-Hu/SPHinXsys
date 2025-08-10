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
 * @file zero_gradient_residual.h
 * @brief The error residual for computing the gradient of unit.
 * @author	Xiangyu Hu
 */

#ifndef ZERO_GRADIENT_RESIDUAL_H
#define ZERO_GRADIENT_RESIDUAL_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class BaseInteractionType>
class ZeroGradientResidualBase : public BaseInteractionType
{
  public:
    template <class DynamicsIdentifier>
    explicit ZeroGradientResidualBase(DynamicsIdentifier &identifier);
    virtual ~ZeroGradientResidualBase() {}

  protected:
    DiscreteVariable<Vecd> *dv_zero_gradient_residue_; ///< "ZeroGradientResidue"
};

template <typename...>
class ZeroGradientResidual;

template <class KernelCorrectionType, typename... Parameters>
class ZeroGradientResidual<Inner<KernelCorrectionType, Parameters...>>
    : public ZeroGradientResidualBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = ZeroGradientResidualBase<Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit ZeroGradientResidual(Inner<Parameters...> &inner_relation);
    virtual ~ZeroGradientResidual() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Vecd *zero_gradient_residue_;
        Real *Vol_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
};

template <class KernelCorrectionType, typename... Parameters>
class ZeroGradientResidual<Contact<Boundary, KernelCorrectionType, Parameters...>>
    : public ZeroGradientResidualBase<Interaction<Contact<Parameters...>>>
{
    using BaseInteraction = ZeroGradientResidualBase<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit ZeroGradientResidual(Contact<Parameters...> &contact_relation);
    virtual ~ZeroGradientResidual() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Vecd *zero_gradient_residue_;
        Real *contact_Vol_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
};
} // namespace SPH
#endif // ZERO_GRADIENT_RESIDUAL_H
