#ifndef SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP
#define SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_variable_array.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, template <typename> class VariableType>
VariableDataArray VariableArray<VariableType<DataType>>::
    DelegatedFieldArray(const ParallelDevicePolicy &par_device)
{
    if (!isFieldArrayDelegated())
    {
        device_shared_variable_array_keeper_
            .createPtr<VariableArray<DeviceShared, VariableType<DataType>>>(this);
    }
    return VariableDataArray(delegated_data_array_, getArraySize());
};
//=================================================================================================//
template <typename DataType, template <typename> class VariableType>
VariableArray<DeviceShared, VariableType<DataType>>::
    VariableArray(VariableType<DataType> *host_variable_array)
    : Entity(host_variable_array->Name()), device_shared_data_array_(nullptr)
{
    StdVec<VariableType<DataType> *> host_variables = host_variable_array->getVariables();
    device_shared_data_array_ = allocateDeviceShared<VariableField<DataType>>(host_variables.size());

    for (size_t i = 0; i != host_variables.size(); ++i)
    {
        device_shared_data_array_[i] = host_variables[i]->DelegatedData(ParallelDevicePolicy{});
    }
    host_variable_array->setDelegateDataArray(device_shared_data_array_);
}
//=================================================================================================//
template <typename DataType>
VariableArray<DeviceShared, VariableType<DataType>>::~VariableArray()
{
    freeDeviceData(device_shared_data_array_);
}
//=================================================================================================//
} // namespace SPH

#endif // SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP