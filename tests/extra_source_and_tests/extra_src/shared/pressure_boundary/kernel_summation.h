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
 * @file 	kernel_summation.h
 * @brief   Here, according to the zeroth order consistency, we calculate the
            kernel summation for imposing the pressure boundary condition.
 * @author	Shuoguo Zhang and Xiangyu Hu
 */

#ifndef KERNEL_SUMMATION_H
#define KERNEL_SUMMATION_H

#include "base_general_dynamics.h"

namespace SPH
{
template <typename... InteractionTypes>
class NablaWV;

template <class DataDelegationType>
class NablaWV<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit NablaWV(BaseRelationType &base_relation);
    virtual ~NablaWV() {};

  protected:
    Vecd *kernel_sum_;
};

template <>
class NablaWV<Inner<>>
    : public NablaWV<DataDelegateInner>
{
  public:
    explicit NablaWV(BaseInnerRelation &inner_relation);
    virtual ~NablaWV() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *Vol_;
};

template <>
class NablaWV<Contact<>>
    : public NablaWV<DataDelegateContact>
{
  public:
    explicit NablaWV(BaseContactRelation &contact_relation)
        : NablaWV<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };
    virtual ~NablaWV() {};
    void interaction(size_t index_i, Real dt = 0.0);

    StdVec<Real *> contact_Vol_;
};

using NablaWVComplex = ComplexInteraction<NablaWV<Inner<>, Contact<>>>;

} // namespace SPH
#endif // KERNEL_SUMMATION_H
