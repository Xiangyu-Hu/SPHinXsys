#ifndef SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP
#define SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_variable_array.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, template <typename> class VariableType>
DataArray<DataType> *VariableArray<DataType, VariableType>::
    DelegatedDataArray(const ParallelDevicePolicy &par_device)
{
    if (!isDataArrayDelegated())
    {
        device_only_variable_array_keeper_
            .createPtr<DeviceOnlyVariableArray<DataType, VariableType>>(this);
    }
    return delegated_data_array_;
}
//=================================================================================================//
template <typename DataType, template <typename> class VariableType>
DeviceOnlyVariableArray<DataType, VariableType>::
    DeviceOnlyVariableArray(VariableArray<DataType, VariableType> *host_variable_array)
    : Entity(host_variable_array->Name()), device_only_data_array_(nullptr)
{
    StdVec<VariableType<DataType> *> host_variables = host_variable_array->getVariables();
    size_t data_size = host_variable_array->getArraySize();
    DataArray<DataType> *host_data = host_variable_array->Data();
    device_only_data_array_ = allocateDeviceOnly<DataArray<DataType>>(data_size);
    for (size_t i = 0; i != data_size; ++i)
    {
        DataType *data = host_variables[i]->DelegatedData(ParallelDevicePolicy{});
        copyToDevice(data, device_only_data_array_ + i, 1);
    }
    host_variable_array->setDelegateDataArray(device_only_data_array_);
}
//=================================================================================================//
template <typename DataType, template <typename> class VariableType>
DeviceOnlyVariableArray<DataType, VariableType>::~DeviceOnlyVariableArray()
{
    freeDeviceData(device_only_data_array_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP
