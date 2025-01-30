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
 * @file density_regularization.h
 * @brief Here, we define the algorithm classes for computing
 * the density of a continuum by kernel function summation.
 * @details We are using templates and their explicit or partial specializations
 * to identify variations of the interaction types..
 * @author Xiangyu Hu
 */

#ifndef KERNEL_SUMMATION_CK_H
#define KERNEL_SUMMATION_CK_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"

namespace SPH
{

namespace fluid_dynamics
{
// template <class BaseInteractionType>
// class NablaWVCKBase : public BaseInteractionType
// {

//   public:
//     template <class DynamicsIdentifier>
//     explicit NablaWVCKBase(DynamicsIdentifier &identifier);
//     virtual ~NablaWVCKBase() {};

//   protected:
//     DiscreteVariable<Vecd> *dv_kernel_sum_;
// };
template <typename... RelationTypes>
class NablaWV;

template <template <typename...> class RelationType, typename... Parameters>
class NablaWV<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{

  public:
    template <class DynamicsIdentifier>
    explicit NablaWV(DynamicsIdentifier &identifier);
    virtual ~NablaWV() {};

    class InteractKernel
        : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       NablaWV<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);

      protected:
        Vecd *kernel_sum_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_kernel_sum_;
};

template <class FlowType, typename... Parameters>
class NablaWV<Inner<FlowType, Parameters...>>
    : public NablaWV<Base, Inner<Parameters...>>
{

  public:
    explicit NablaWV(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~NablaWV() {};

    class InteractKernel
        : public NablaWV<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       NablaWV<Inner<FlowType, Parameters...>> &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol_;
    };

  protected:
    DiscreteVariable<Real> *dv_Vol_;
};

template <typename... Parameters>
class NablaWV<Contact<Parameters...>>
    : public NablaWV<Base, Contact<Parameters...>>
{

  public:
    explicit NablaWV(Relation<Contact<Parameters...>> &contact_relation);
    virtual ~NablaWV() {};

    class InteractKernel
        : public NablaWV<Base, Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       NablaWV<Contact<Parameters...>> &encloser,
                       size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *contact_wall_Vol_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_wall_Vol_;
};
//--------------------------------------------------------------------------------------
// Alias Definitions for Specific Configurations
//--------------------------------------------------------------------------------------
using NablaWVInnerCK = NablaWV<Inner<Internal>>;
using NablaWVComplexCK = NablaWV<Inner<Internal>, Contact<>>;
using NablaWVComplexFreeSurfaceCK = NablaWV<Inner<FreeSurface>, Contact<>>;

} // namespace fluid_dynamics
} // namespace SPH
#endif // KERNEL_SUMMATION_CK_H
