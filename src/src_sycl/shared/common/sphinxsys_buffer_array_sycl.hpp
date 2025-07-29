#ifndef SPHINXSYS_BUFFER_ARRAY_SYCL_HPP
#define SPHINXSYS_BUFFER_ARRAY_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_buffer_array.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DataArray<DataType> *BufferArray<DataType>::
    DelegatedDataArray(const DeviceExecution<PolicyType> &ex_policy)
{
    if (!isDataArrayDelegated())
    {
        device_buffer_array_keeper_
            .createPtr<DeviceBufferArray<DataType>>(ex_policy, this);
    }
    return delegated_data_array_;
}
//=================================================================================================//
template <typename DataType>
void BufferArray<DataType>::setHostStagingDataArray()
{
    if (!isDataArrayHostStaged())
    {
        host_staging_data_array_ =
            device_buffer_array_keeper_.getPtr()->allocateHostStagingDataArray();
    }
}
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DataArray<DataType> *BufferArray<DataType>::
    getHostStagingDataArray(const DeviceExecution<PolicyType> &ex_policy, UnsignedInt data_size)
{
    for (size_t i = 0; i != array_size_; ++i)
    {
        copyFromDevice(host_staging_data_array_[i], delegated_data_array_[i], data_size);
    }
    return host_staging_data_array_;
}
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DeviceBufferArray<DataType>::
    DeviceBufferArray(const DeviceExecution<PolicyType> &ex_policy,
                      BufferArray<DataType> *host_buffer_array)
    : Entity(host_buffer_array->Name()), array_size_(host_buffer_array->getArraySize()),
      device_shared_data_array_(allocateDeviceOnly<DataArray<DataType>>(array_size_)),
      host_staging_data_array_(nullptr)
{
    UnsignedInt buffer_size = host_buffer_array->getBufferSize();
    for (size_t i = 0; i != array_size_; ++i)
    {
        device_shared_data_array_[i] = allocateDeviceOnly<DataType>(buffer_size);
    }
    host_buffer_array->setDelegateDataArray(device_shared_data_array_);
}
//=================================================================================================//
template <typename DataType>
DataArray<DataType> *DeviceBufferArray<DataType>::allocateHostStagingDataArray()
{
    if (host_staging_data_array_ == nullptr)
    {
        host_staging_data_array_ = allocateHostStaging<DataArray<DataType>>(array_size_);
        for (size_t i = 0; i != array_size_; ++i)
        {
            host_staging_data_array_[i] = allocateHostStaging<DataType>(getBufferSize());
        }
    }
    return host_staging_data_array_;
}
//=================================================================================================//
template <typename DataType>
DeviceBufferArray<DataType>::~DeviceBufferArray()
{
    for (size_t i = 0; i != array_size_; ++i)
    {
        freeDeviceData(device_shared_data_array_[i]);
        if (host_staging_data_array_[i] != nullptr)
        {
            freeDeviceData(host_staging_data_array_[i]);
        }
    }
    freeDeviceData(device_shared_data_array_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_BUFFER_ARRAY_SYCL_HPP
