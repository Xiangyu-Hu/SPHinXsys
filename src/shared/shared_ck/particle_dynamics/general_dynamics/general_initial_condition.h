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
 * @file general_initial_condition.h
 * @brief tbd.
 * @author Xiangyu Hu
 */

#ifndef GENERAL_INITIAL_CONDITION_H
#define GENERAL_INITIAL_CONDITION_H

#include "base_general_dynamics.h"

namespace SPH
{
template <class DynamicsIdentifier, typename InitializeFunctionType>
class InitialCondition : public BaseLocalDynamics<DynamicsIdentifier>
{
    using DataType = typename InitializeFunctionType::ReturnType;
    using Initialize = typename InitializeFunctionType::ComputingKernel;

  public:
    template <typename... Args>
    InitialCondition(DynamicsIdentifier &identifier, const std::string &variable_name, Args &...args)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
          dv_variable_(this->particles_->template registerStateVariableOnly<DataType>(variable_name)),
          initialize_method_(this->particles_, std::forward<Args>(args)...){};
    virtual ~InitialCondition() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
              initialize_(ex_policy, encloser.initialize_method_){};

        void update(size_t index_i, Real dt = 0.0)
        {
            variable_[index_i] = initialize_(pos_[index_i]);
        };

      protected:
        Vecd *pos_;
        DataType *variable_;
        Initialize initialize_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<DataType> *dv_variable_;
    InitializeFunctionType initialize_method_;
};

template <typename DataType>
class UniformDistribution : public ReturnFunction<DataType>
{
  public:
    UniformDistribution(BaseParticles *particles, DataType constant_value)
        : constant_value_(constant_value) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : constant_value_(encloser.constant_value_){};

        DataType operator()(const Vecd &position) { return constant_value_; };

      protected:
        DataType constant_value_;
    };

  protected:
    DataType constant_value_;
};
} // namespace SPH
#endif // GENERAL_INITIAL_CONDITION_H
