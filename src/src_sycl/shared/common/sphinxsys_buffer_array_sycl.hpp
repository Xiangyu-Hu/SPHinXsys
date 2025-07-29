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
        device_shared_buffer_array_keeper_
            .createPtr<DeviceSharedBufferArray<DataType>>(ex_policy, this);
    }
    return delegated_data_array_;
}
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DeviceSharedBufferArray<DataType>::
    DeviceSharedBufferArray(const DeviceExecution<PolicyType> &ex_policy,
                            BufferArray<DataType> *host_buffer_array)
    : Entity(host_buffer_array->Name()), array_size_(host_buffer_array->getArraySize()),
      device_shared_data_array_(allocateDeviceShared<DataArray<DataType>>(array_size_))
{
    UnsignedInt buffer_size = host_buffer_array->getBufferSize();
    for (size_t i = 0; i != array_size_; ++i)
    {
        device_shared_data_array_[i] = allocateDeviceShared<DataType>(buffer_size);
    }
    host_buffer_array->setDelegateDataArray(device_shared_data_array_);
}
//=================================================================================================//
template <typename DataType>
DeviceSharedBufferArray<DataType>::~DeviceSharedBufferArray()
{
    for (size_t i = 0; i != array_size_; ++i)
    {
        freeDeviceData(device_shared_data_array_[i]);
    }
    freeDeviceData(device_shared_data_array_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_BUFFER_ARRAY_SYCL_HPP
