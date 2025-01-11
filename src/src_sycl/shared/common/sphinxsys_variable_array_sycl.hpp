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
    device_only_data_array_ = allocateDeviceOnly<DataArray<DataType>>(host_variables.size());

    for (size_t i = 0; i != host_variables.size(); ++i)
    {
        device_only_data_array_[i] = host_variables[i]->DelegatedData(ParallelDevicePolicy{});
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
