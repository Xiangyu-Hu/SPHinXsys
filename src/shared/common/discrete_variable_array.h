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
 * @file discrete_variable_array.h
 * @brief tbd
 * @author Xiangyu Hu
 */

#ifndef DISCRETE_VARIABLE_ARRAY_H
#define DISCRETE_VARIABLE_ARRAY_H

#include "sphinxsys_variable.h"

namespace SPH
{

template <typename DataType>
using VariableField = DataType *;

template <typename DataType>
class DiscreteVariableArray : public Entity
{
    UniquePtrKeeper<Entity> device_shared_discrete_variable_array_keeper_;

  public:
    DiscreteVariableArray(const std::string &name, StdVec<DiscreteVariable<DataType> *> discrete_variables)
        : Entity(name), discrete_variables_(discrete_variables),
          field_array_(nullptr), delegated_field_array_(nullptr)
    {
        field_array_ = new VariableField<DataType>[discrete_variables.size()];
        for (size_t i = 0; i != discrete_variables.size(); ++i)
        {
            field_array_[i] = discrete_variables[i]->DataField();
        }
        delegated_field_array_ = field_array_;
    };
    ~DiscreteVariableArray() { delete[] field_array_; };
    StdVec<DiscreteVariable<DataType> *> getDiscreteVariables() { return discrete_variables_; };
    DataType *VariableArray() { return field_array_; };

    template <class ExecutionPolicy>
    VariableField<DataType> *DelegatedFieldArray(const ExecutionPolicy &ex_policy) { return field_array_; };
    VariableField<DataType> *DelegatedFieldArray(const ParallelDevicePolicy &par_device);
    size_t getFieldArraySize() { return discrete_variables_.size(); }

    bool isFieldArrayDelegated() { return field_array_ != delegated_field_array_; };
    void setDelegateFieldArray(VariableField<DataType> *field_array_) { delegated_field_array_ = field_array_; };

  private:
    StdVec<DiscreteVariable<DataType> *> discrete_variables_;
    VariableField<DataType> *field_array_;
    VariableField<DataType> *delegated_field_array_;
};

template <typename DataType>
class DeviceSharedDiscreteVariableArray : public Entity
{
  public:
    DeviceSharedDiscreteVariableArray(DiscreteVariableArray<DataType> *host_variable_array);
    ~DeviceSharedDiscreteVariableArray();

  protected:
    VariableField<DataType> *device_shared_field_array_;
};
} // namespace SPH
#endif // DISCRETE_VARIABLE_ARRAY_H
