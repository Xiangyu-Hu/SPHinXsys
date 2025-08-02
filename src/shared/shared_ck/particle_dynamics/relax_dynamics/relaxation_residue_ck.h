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
 * @file relaxation_residue_ck.h
 * @brief TBD.
 * @author Xiangyu Hu
 */

#ifndef RELAXATION_RESIDUE_CK_H
#define RELAXATION_RESIDUE_CK_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class BaseInteractionType>
class RelaxationResidueBase : public BaseInteractionType
{
  public:
    template <class DynamicsIdentifier>
    explicit RelaxationResidueBase(DynamicsIdentifier &identifier);
    virtual ~RelaxationResidueBase() {}

  protected:
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Vecd> *dv_residue_;
};

template <typename...>
class RelaxationResidueCK;

template <class KernelCorrectionType, typename... Parameters>
class RelaxationResidueCK<Inner<KernelCorrectionType, Parameters...>>
    : public RelaxationResidueBase<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = RelaxationResidueBase<Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit RelaxationResidueCK(Inner<Parameters...> &inner_relation);
    virtual ~RelaxationResidueCK() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Real *Vol_;
        Vecd *residue_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
};

template <class KernelCorrectionType, typename... Parameters>
class RelaxationResidueCK<Contact<Boundary, KernelCorrectionType, Parameters...>>
    : public RelaxationResidueBase<Interaction<Contact<Parameters...>>>
{
    using BaseInteraction = RelaxationResidueBase<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit RelaxationResidueCK(Contact<Parameters...> &contact_relation);
    virtual ~RelaxationResidueCK() {}

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, 
          EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        Real *contact_Vol_;
        Vecd *residue_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
};
} // namespace SPH
#endif // RELAXATION_RESIDUE_CK_H
