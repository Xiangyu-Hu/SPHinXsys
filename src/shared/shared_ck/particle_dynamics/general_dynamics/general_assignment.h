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
 * @file general_assignment.h
 * @brief tbd.
 * @author Xiangyu Hu
 */

#ifndef GENERAL_ASSIGNMENT_H
#define GENERAL_ASSIGNMENT_H

#include "base_local_dynamics.h"

#include <string>
#include <utility>

namespace SPH
{
template <class DynamicsIdentifier, typename AssignmentFunctionType>
class VariableAssignment : public BaseLocalDynamics<DynamicsIdentifier>
{
    using Assign = typename AssignmentFunctionType::ComputingKernel;

  public:
    template <typename... Args>
    VariableAssignment(DynamicsIdentifier &identifier, Args &&...args);
    virtual ~VariableAssignment() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0) { assign_(index_i); };

      protected:
        Assign assign_;
    };

  protected:
    AssignmentFunctionType assignment_method_;
};

template <typename...>
class SpatialDistribution;

template <typename DistributionType>
class SpatialDistribution<DistributionType>
    : public ReturnFunction<typename DistributionType::ReturnType>
{
    using DataType = typename DistributionType::ReturnType;

  public:
    template <typename... Args>
    SpatialDistribution(BaseParticles *particles, const std::string &variable_name, Args &&...args);

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void operator()(UnsignedInt index_i) { variable_[index_i] = distribution_(pos_[index_i]); };

      protected:
        DataView<DataType> variable_;
        DataView<Vecd> pos_;
        DistributionType distribution_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
    DiscreteVariable<Vecd> *dv_pos_;
    DistributionType distribution_;
};

template <typename DataType>
class ConstantValue
{
  public:
    ConstantValue(BaseParticles *particles, const std::string &variable_name, DataType constant_value);

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void operator()(UnsignedInt index_i) { variable_[index_i] = constant_value_; };

      protected:
        DataView<DataType> variable_;
        DataType constant_value_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
    DataType constant_value_;
};
} // namespace SPH
#endif // GENERAL_ASSIGNMENT_H
