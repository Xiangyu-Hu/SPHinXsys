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
 * @file general_assignment.h
 * @brief tbd.
 * @author Xiangyu Hu
 */

#ifndef GENERAL_ASSIGNMENT_H
#define GENERAL_ASSIGNMENT_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class DynamicsIdentifier, typename AssignmentFunctionType>
class VariableAssignment : public BaseLocalDynamics<DynamicsIdentifier>
{
    using DataType = typename AssignmentFunctionType::ReturnType;
    using Assign = typename AssignmentFunctionType::ComputingKernel;

  public:
    template <typename... Args>
    VariableAssignment(DynamicsIdentifier &identifier, const std::string &variable_name, Args &&...args)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          dv_variable_(this->particles_->template registerStateVariableOnly<DataType>(variable_name)),
          assignment_method_(this->particles_, std::forward<Args>(args)...){};
    virtual ~VariableAssignment() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
              assign_(ex_policy, encloser.assignment_method_){};

        void update(UnsignedInt index_i, Real dt = 0.0)
        {
            variable_[index_i] = assign_(index_i);
        };

      protected:
        DataType *variable_;
        Assign assign_;
    };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
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
    SpatialDistribution(BaseParticles *particles, Args &&...args)
        : dv_pos_(particles->template getVariableByName<Vecd>("Position")),
          distribution_(std::forward<Args>(args)...){};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              distribution_(encloser.distribution_){};

        DataType operator()(UnsignedInt index_i) { return distribution_(pos_[index_i]); };

      protected:
        Vecd *pos_;
        DistributionType distribution_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    DistributionType distribution_;
};

template <typename DataType>
class ConstantValue : public ReturnFunction<DataType>
{
  public:
    ConstantValue(BaseParticles *particles, DataType constant_value)
        : constant_value_(constant_value) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : constant_value_(encloser.constant_value_){};

        DataType operator()(UnsignedInt index_i) { return constant_value_; };

      protected:
        DataType constant_value_;
    };

  protected:
    DataType constant_value_;
};

template <typename DataType>
class CopyVariable : public ReturnFunction<DataType>
{
  public:
    CopyVariable(BaseParticles *particles, const std::string &copy_variable_name)
        : dv_copy_variable_(particles->template getVariableByName<DataType>(
              copy_variable_name)) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : copy_variable_(encloser.dv_copy_variable_->DelegatedData(ex_policy)){};

        DataType operator()(UnsignedInt index_i) { return copy_variable_[index_i]; };

      protected:
        DataType *copy_variable_;
    };

  protected:
    DiscreteVariable<DataType> *dv_copy_variable_;
};
} // namespace SPH
#endif // GENERAL_ASSIGNMENT_H
