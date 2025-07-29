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
 * @file sphinxsys_buffer_array.h
 * @brief tbd
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_BUFFER_ARRAY_H
#define SPHINXSYS_BUFFER_ARRAY_H

#include "sphinxsys_variable_array.h"

namespace SPH
{
template <typename DataType>
class BufferArray : public Entity
{
    UniquePtrKeeper<Entity> device_shared_buffer_array_keeper_;

  public:
    BufferArray(StdVec<DiscreteVariable<DataType> *> variables, UnsignedInt buffer_size)
        : Entity("BufferArray"), variables_(variables),
          buffer_size_(buffer_size), data_array_(nullptr), delegated_data_array_(nullptr)
    {
        data_array_ = new DataArray<DataType>[variables.size()];
        for (size_t i = 0; i != variables.size(); ++i)
        {
            data_array_[i] = new DataType[buffer_size];
        }
        delegated_data_array_ = data_array_;
    };
    ~BufferArray() { delete[] data_array_; };
    StdVec<DiscreteVariable<DataType> *> getVariables() { return variables_; };
    size_t getArraySize() { return variables_.size(); }
    UnsignedInt getBufferSize() { return buffer_size_; }
    DataArray<DataType> *Data() { return data_array_; };
    bool isDataArrayDelegated() { return data_array_ != delegated_data_array_; };

    template <class ExecutionPolicy>
    DataArray<DataType> *DelegatedDataArray(const ExecutionPolicy &ex_policy)
    {
        return data_array_;
    };

    template <class PolicyType>
    DataArray<DataType> *DelegatedDataArray(const DeviceExecution<PolicyType> &ex_policy);

    void setDelegateDataArray(DataArray<DataType> *data_array_)
    {
        delegated_data_array_ = data_array_;
    };

  private:
    StdVec<DiscreteVariable<DataType> *> variables_;
    UnsignedInt buffer_size_;
    DataArray<DataType> *data_array_;
    DataArray<DataType> *delegated_data_array_;
};

template <typename DataType>
class DeviceSharedBufferArray : public Entity
{
  public:
    template <class PolicyType>
    DeviceSharedBufferArray(const DeviceExecution<PolicyType> &ex_policy,
                            BufferArray<DataType> *host_buffer_array);
    ~DeviceSharedBufferArray();

  protected:
    size_t array_size_;
    DataArray<DataType> *device_shared_data_array_;
};
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_H
