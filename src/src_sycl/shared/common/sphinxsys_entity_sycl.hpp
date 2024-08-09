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
 * @file    sphinxsys_entity_sycl.hpp
 * @brief   TBD.
 * @author  Alberto Guarnieri and Xiangyu Hu
 */
#ifndef SPHINXSYS_ENTITY_SYCL_HPP
#define SPHINXSYS_ENTITY_SYCL_HPP

#include "execution_sycl.h"
#include "sphinxsys_entity.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DeviceSharedSingularVariable<DataType>::
    DeviceSharedSingularVariable(SingularVariable<DataType> *host_variable)
    : BaseEntity(host_variable->Name()),
      device_shared_value_(allocateDeviceShared<DataType>(1))
{
    *device_shared_value_ = *host_variable->ValueAddress();
    host_variable->setDelegateValueAddress(device_shared_value_);
}
//=================================================================================================//
template <typename DataType>
DeviceSharedSingularVariable<DataType>::~DeviceSharedSingularVariable()
{
    freeDeviceData(device_shared_value_);
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyConstantEntity<DataType>::DeviceOnlyConstantEntity(ConstantEntity<DataType> *host_constant)
    : BaseEntity(host_constant->Name()), device_only_value_(allocateDeviceOnly<DataType>(1))
{
    *device_only_value_ = *host_constant->ValueAddress();
    host_constant->setDelegateValueAddress(device_only_value_);
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyConstantEntity<DataType>::~DeviceOnlyConstantEntity()
{
    freeDeviceData(device_only_value_);
}
//=================================================================================================//
template <typename DataType>
void DiscreteVariable<DataType>::synchronizeWithDevice()
{
    if (existDeviceDataField())
    {
        copyFromDevice(data_field_, device_data_field_, data_size_);
    }
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyDiscreteVariable<DataType>::
    DeviceOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable)
    : BaseEntity(host_variable->Name()), device_only_data_field_(nullptr)
{
    size_t data_size = host_variable->getDataSize();
    device_only_data_field_ = allocateDeviceOnly<DataType>(data_size);
    copyToDevice(host_variable->DataField(), device_only_data_field_, data_size);
    host_variable->setDeviceDataField(device_only_data_field_);
};
//=================================================================================================//
template <typename DataType>
DeviceOnlyDiscreteVariable<DataType>::~DeviceOnlyDiscreteVariable()
{
    freeDeviceData(device_only_data_field_);
}
//=================================================================================================//
} // namespace SPH

#endif // SPHINXSYS_ENTITY_SYCL_HPP