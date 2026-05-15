#ifndef SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP
#define SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_variable_array.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DataPtr<DataType> *VariableArray<DataType>::DelegatedOnDevice()
{
    if (!isDataArrayDelegated())
    {
        device_only_variable_array_ = device_only_variable_array_keeper_.createPtr<
            DeviceOnlyVariableArray<DataType>>(DeviceExecution<PolicyType>{}, this);
    }
    return device_only_variable_array_->DeviceOnlyDataPtr();
}
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DeviceOnlyVariableArray<DataType>::
    DeviceOnlyVariableArray(const DeviceExecution<PolicyType> &ex_policy,
                            VariableArray<DataType> *host_variable_array)
    : Quantity(host_variable_array->Name()), device_only_data_ptr_(nullptr)
{
    StdVec<DiscreteVariable<DataType> *> host_variables = host_variable_array->getVariables();
    size_t data_size = host_variable_array->getArraySize();
    device_only_data_ptr_ = allocateDeviceOnly<DataPtr<DataType>>(data_size);
    for (size_t i = 0; i != data_size; ++i)
    {
        DataType *data = host_variables[i]->DelegatedData(ex_policy);
        copyToDevice(data, device_only_data_ptr_ + i, 1);
    }
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyVariableArray<DataType>::~DeviceOnlyVariableArray()
{
    freeDeviceData(device_only_data_ptr_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP
