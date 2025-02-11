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
 * @file 	interpolation_dynamics.h
 * @brief 	There are the classes for interpolation algorithm.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef INTERPOLATION_DYNAMICS_H
#define INTERPOLATION_DYNAMICS_H

#include "interaction_algorithms_ck.hpp"

namespace SPH
{

template <typename...>
class Interpolation;
/**
 * @class Interpolation
 * @brief Base class for interpolation.
 */
template <typename DataType>
class Interpolation<Contact<DataType>> : public Interaction<Contact<>>
{
  public:
    Interpolation(Relation<Contact<>> &pair_contact_relation, const std::string &variable_name);
    virtual ~Interpolation() {};

    class InteractKernel : public Interaction<Contact<>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        DataType zero_value_;
        DataType *interpolated_quantities_;
        Real *contact_Vol_;
        DataType *contact_data_;
    };

  protected:
    DiscreteVariable<DataType> *dv_interpolated_quantities_;
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
    StdVec<DiscreteVariable<DataType> *> dv_contact_data_;
};

template <class ExecutionPolicy, typename DataType>
class ObservingAQuantityCK : public InteractionDynamicsCK<ExecutionPolicy, Interpolation<Contact<DataType>>>
{
  public:
    ObservingAQuantityCK(Relation<Contact<>> &pair_contact_relation, const std::string &variable_name)
        : InteractionDynamicsCK<ExecutionPolicy, Interpolation<Contact<DataType>>>(pair_contact_relation, variable_name) {};
    template <typename BodyRelationType, typename FirstArg>
    explicit ObservingAQuantityCK(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : ObservingAQuantityCK(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~ObservingAQuantityCK() {};
};
} // namespace SPH
#endif // INTERPOLATION_DYNAMICS_H
