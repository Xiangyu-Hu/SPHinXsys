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
 * @brief A buffer array defines an device shared array of temporary data corresponding to
 * a discrete variable array so that the data in the array can be temporarily saved
 * (at both device or host) in the buffer. I expect that these temporary data
 * can be manipulated by both host and device as they are device shared.
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_BUFFER_ARRAY_H
#define SPHINXSYS_BUFFER_ARRAY_H

#include "sphinxsys_variable_array.h"

namespace SPH
{

class DeviceBufferArray;

template <typename DataType>
class BufferArray : public Entity
{
    UniquePtrKeeper<DeviceBufferArray> device_buffer_array_keeper_;

  public:
    BufferArray(StdVec<DiscreteVariable<DataType> *> variables, UnsignedInt buffer_size)
        : Entity("BufferArray"), variables_(variables),
          buffer_size_(buffer_size), data_array_(nullptr),
          delegated_data_array_(nullptr), host_staging_data_array_(nullptr)
    {
        data_array_ = new DataArray<DataType>[variables.size()];
        for (size_t i = 0; i != variables.size(); ++i)
        {
            data_array_[i] = new DataType[buffer_size];
        }
        delegated_data_array_ = data_array_;
        host_staging_data_array_ = data_array_;
    };

    ~BufferArray()
    {
        for (size_t i = 0; i != getArraySize(); ++i)
        {
            delete[] data_array_[i];
        }
        delete[] data_array_;
    };

    StdVec<DiscreteVariable<DataType> *> getVariables() { return variables_; };
    size_t getArraySize() { return variables_.size(); }
    UnsignedInt getBufferSize() { return buffer_size_; }
    DataArray<DataType> *Data() { return data_array_; };
    bool isDataArrayDelegated() { return data_array_ != delegated_data_array_; };
    bool isDataArrayHostStaged() { return data_array_ != host_staging_data_array_; };

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

    void setHostStagingDataArray();

    template <class ExecutionPolicy>
    DataArray<DataType> *getHostStagingDataArray(
        const ExecutionPolicy &ex_policy, UnsignedInt data_size)
    {
        return host_staging_data_array_;
    };

    template <class PolicyType>
    DataArray<DataType> *getHostStagingDataArray(
        const DeviceExecution<PolicyType> &ex_policy, UnsignedInt data_size);

  private:
    StdVec<DiscreteVariable<DataType> *> variables_;
    UnsignedInt buffer_size_;
    DataArray<DataType> *data_array_;
    DataArray<DataType> *delegated_data_array_;
    DataArray<DataType> *host_staging_data_array_;
};

template <typename DataType>
class DeviceBufferArray : public Entity
{
  public:
    template <class PolicyType>
    DeviceBufferArray(const DeviceExecution<PolicyType> &ex_policy,
                      BufferArray<DataType> *host_buffer_array);
    ~DeviceBufferArray();

    DataArray<DataType> *allocateHostStagingDataArray();

  protected:
    size_t array_size_;
    DataArray<DataType> *device_only_data_array_;
    DataArray<DataType> *host_staging_data_array_;
};

template <typename DataType>
using AllocatedDataArrayPair = std::pair<AllocatedDataArray<DataType>, AllocatedDataArray<DataType>>;

template <typename DataType>
using AllocatedDataArrayPairSet = std::pair<AllocatedDataArrayPair<DataType>, UnsignedInt>;

struct CopyAllocatedDataArrayPairSet
{
    template <typename DataType>
    void operator()(AllocatedDataArrayPairSet<DataType> &allocated_pair_set,
                    size_t source_data_index, size_t target_data_index)
    {
        auto &source_allocation = allocated_pair_set.first.first;
        auto &target_allocation = allocated_pair_set.first.second;
        for (size_t i = 0; i < allocated_pair_set.second; ++i)
        {
            target_allocation[i][target_data_index] = source_allocation[i][source_data_index];
        }
    }
};
} // namespace SPH
#endif // SPHINXSYS_BUFFER_ARRAY_H
