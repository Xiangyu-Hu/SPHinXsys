#ifndef DISCRETE_VARIABLE_ARRAY_SYCL_HPP
#define DISCRETE_VARIABLE_ARRAY_SYCL_HPP

#include "discrete_variable_array.h"
#include "implementation_sycl.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
VariableField<DataType> *DiscreteVariableArray<DataType>::
    DelegatedFieldArray(const ParallelDevicePolicy &par_device)
{
    if (!isFieldArrayDelegated())
    {
        device_shared_discrete_variable_array_keeper_
            .createPtr<DeviceSharedDiscreteVariableArray<DataType>>(this);
    }
    return delegated_field_array_;
};
//=================================================================================================//
template <typename DataType>
DeviceSharedDiscreteVariableArray<DataType>::
    DeviceSharedDiscreteVariableArray(DiscreteVariableArray<DataType> *host_variable_array)
    : Entity(host_variable_array->Name()), device_shared_field_array_(nullptr)
{
    StdVec<DiscreteVariable<DataType> *> host_discrete_variables = host_variable_array->getDiscreteVariables();
    device_shared_field_array_ = allocateDeviceShared<VariableField<DataType>>(host_discrete_variables.size());

    for (size_t i = 0; i != host_discrete_variables.size(); ++i)
    {
        device_shared_field_array_[i] = host_discrete_variables[i]->DelegatedData(ParallelDevicePolicy{});
    }
    host_variable_array->setDelegateValueAddress(device_shared_field_array_);
}
//=================================================================================================//
template <typename DataType>
DeviceSharedSingularVariable<DataType>::~DeviceSharedSingularVariable()
{
    freeDeviceData(device_shared_field_array_);
}
//=================================================================================================//
} // namespace SPH

#endif // DISCRETE_VARIABLE_ARRAY_SYCL_HPP