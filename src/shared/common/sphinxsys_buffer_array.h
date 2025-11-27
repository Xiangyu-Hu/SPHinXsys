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
 * @file sphinxsys_buffer_array.h
 * @brief A buffer array defines temporary data corresponding to
 * a discrete variable array so that the data in the array can be temporarily saved
 * at device in the buffer.
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_BUFFER_ARRAY_H
#define SPHINXSYS_BUFFER_ARRAY_H

#include "sphinxsys_variable_array.h"

namespace SPH
{
template <typename DataType>
class DeviceVariableBufferArray;

template <typename DataType>
class VariableBufferArray : public DiscreteVariableArray<DataType>
{
    UniquePtrKeeper<DeviceVariableBufferArray<DataType>> device_buffer_array_keeper_;

  public:
    VariableBufferArray(StdVec<DiscreteVariable<DataType> *> variables, UnsignedInt buffer_size)
        : DiscreteVariableArray<DataType>(variables), buffer_size_(buffer_size),
          buffer_data_ptr_(new DataType[buffer_size * this->array_size_]),
          buffer_array_(new DataArray<DataType>[this->array_size_]),
          delegated_buffer_array_(nullptr), host_staging_buffer_array_(nullptr)
    {
        for (size_t i = 0; i != this->array_size_; ++i)
        {
            buffer_array_[i] = buffer_data_ptr_ + i * buffer_size_;
        }
        delegated_buffer_array_ = buffer_array_;
        host_staging_buffer_array_ = buffer_array_;
    };

    ~VariableBufferArray()
    {
        delete[] buffer_data_ptr_;
        delete[] buffer_array_;
    };

    UnsignedInt getBufferSize() { return buffer_size_; }
    DataArray<DataType> *BufferData() { return buffer_array_; };
    bool isBufferArrayDelegated() { return buffer_array_ != delegated_buffer_array_; };
    bool isBufferArrayHostStaged() { return buffer_array_ != host_staging_buffer_array_; };

    template <class ExecutionPolicy>
    DataArray<DataType> *DelegatedBufferArray(const ExecutionPolicy &ex_policy)
    {
        return buffer_array_;
    };
    template <class PolicyType>
    DataArray<DataType> *DelegatedBufferArrayOnDevice();
    template <class PolicyType>
    DataArray<DataType> *DelegatedBufferArray(const DeviceExecution<PolicyType> &ex_policy)
    {
        return DelegatedBufferArrayOnDevice<PolicyType>();
    };

    void setDelegateBufferArray(DataArray<DataType> *buffer_array_)
    {
        delegated_buffer_array_ = buffer_array_;
    };

    template <class ExecutionPolicy>
    DataArray<DataType> *synchronizeHostStagingBufferArray(
        const ExecutionPolicy &ex_policy, UnsignedInt data_size)
    {
        return host_staging_buffer_array_;
    };

    DataArray<DataType> *synchronizeBufferArrayFromDevice(UnsignedInt data_size);
    template <class PolicyType>
    DataArray<DataType> *synchronizeHostStagingBufferArray(
        const DeviceExecution<PolicyType> &ex_policy, UnsignedInt data_size)
    {
        return synchronizeBufferArrayFromDevice(data_size);
    };

    void setHostStagingBufferArray(DataArray<DataType> *buffer_array_)
    {
        host_staging_buffer_array_ = buffer_array_;
    };

  protected:
    UnsignedInt buffer_size_;
    DataType *buffer_data_ptr_;
    DataArray<DataType> *buffer_array_;
    DataArray<DataType> *delegated_buffer_array_;
    DataArray<DataType> *host_staging_buffer_array_;

  public:
    class CopyVariableToBuffer
    {
        UnsignedInt array_size_;
        DataArray<DataType> *delegated_data_array_;
        DataArray<DataType> *delegated_buffer_array_;

      public:
        template <class ExecutionPolicy>
        CopyVariableToBuffer(const ExecutionPolicy &ex_policy, VariableBufferArray &variable_buffer_array)
            : array_size_(variable_buffer_array.array_size_),
              delegated_data_array_(variable_buffer_array.DelegatedDataArray(ex_policy)),
              delegated_buffer_array_(variable_buffer_array.DelegatedBufferArray(ex_policy)) {}
        void operator()(UnsignedInt source_data_index, UnsignedInt target_data_index)
        {
            for (size_t i = 0; i < array_size_; ++i)
            {
                delegated_buffer_array_[i][target_data_index] =
                    delegated_data_array_[i][source_data_index];
            }
        }
    };
};

template <typename DataType>
class DeviceVariableBufferArray : public Entity
{
  public:
    template <class PolicyType>
    DeviceVariableBufferArray(const DeviceExecution<PolicyType> &ex_policy,
                              VariableBufferArray<DataType> *host_buffer_array);
    ~DeviceVariableBufferArray();
    void allocateHostStagingBufferArray(VariableBufferArray<DataType> *host_buffer_array);
    void synchronizeHostStagingBufferArray(UnsignedInt data_size);

  protected:
    size_t array_size_;
    DataType *device_only_ptr_;
    DataType *host_staging_ptr_;
    UnsignedInt buffer_size_;
    DataArray<DataType> *device_only_buffer_array_;
    DataArray<DataType> *host_staging_buffer_array_;
};
} // namespace SPH
#endif // SPHINXSYS_BUFFER_ARRAY_H
