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
 * @file sphinxsys_variable_array.h
 * @brief tbd
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_VARIABLE_ARRAY_H
#define SPHINXSYS_VARIABLE_ARRAY_H

#include "sphinxsys_variable.h"

namespace SPH
{
template <typename DataType>
using VariableData = DataType *;

template <typename DataType>
using VariableDataArray = std::pair<VariableData *, UnsignedInt>;

class DeviceShared;

template <typename...>
class VariableArray;

template <typename DataType, template <typename> class VariableType>
class VariableArray<DeviceShared, VariableType<DataType>> : public Entity
{
  public:
    VariableArray(VariableArray<VariableType<DataType>> *host_variable_array);
    ~VariableArray();

  protected:
    VariableData<DataType> *device_shared_data_array_;
};

template <typename DataType, template <typename> class VariableType>
class VariableArray<VariableType<DataType>> : public Entity
{
    UniquePtrKeeper<Entity> device_shared_variable_array_keeper_;

  public:
    VariableArray(const std::string &name, StdVec<VariableType<DataType> *> variables)
        : Entity(name), variables_(variables),
          data_array_(nullptr), delegated_data_array_(nullptr)
    {
        data_array_ = new VariableData<DataType>[variables.size()];
        for (size_t i = 0; i != variables.size(); ++i)
        {
            data_array_[i] = variables[i]->Data();
        }
        delegated_data_array_ = data_array_;
    };
    ~DiscreteVariableArray() { delete[] data_array_; };
    StdVec<VariableType<DataType> *> getVariables() { return variables_; };
    DataType *VariableArray() { return data_array_; };

    template <class ExecutionPolicy>
    VariableDataArray DelegatedVariableDataArray(const ExecutionPolicy &ex_policy)
    {
        return VariableDataArray(data_array_, getArraySize());
    };
    VariableDataArray DelegatedVariableDataArray(const ParallelDevicePolicy &par_device);
    size_t getArraySize() { return variables_.size(); }

    bool isFieldArrayDelegated() { return data_array_ != delegated_data_array_; };
    void setDelegateDataArray(VariableData<DataType> *data_array_)
    {
        delegated_data_array_ = data_array_;
    };

  private:
    StdVec<VariableType<DataType> *> variables_;
    VariableData<DataType> *data_array_;
    VariableData<DataType> *delegated_data_array_;
};
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_H
