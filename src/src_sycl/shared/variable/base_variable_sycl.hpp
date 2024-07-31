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
 * @file    base_variable_sycl.hpp
 * @brief   TBD.
 * @author  Alberto Guarnieri and Xiangyu Hu
 */
#ifndef BASE_VARIABLE_SYCL_HPP
#define BASE_VARIABLE_SYCL_HPP

#include "base_variable_sycl.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
SingularDeviceSharedVariable<DataType>::
    SingularVariable(SingularVariable<DataType> &original_variable)
    : BaseVariable(name), device_shared_value_(allocateDeviceShared<DataType>(1))
{
    *device_shared_value_ = *original_variable.ValueAddress();
    original_variable.setDelegateValue(device_shared_value_);
}
//=================================================================================================//
template <typename DataType>
SingularDeviceSharedVariable<DataType>::~SingularDeviceSharedVariable()
{
    freeDeviceData(device_shared_value_);
}
//=================================================================================================//
template <typename DataType>
DiscreteDeviceSharedVariable<DataType>::
    DiscreteVariable(DiscreteVariable<DataType> &original_variable)
    : BaseVariable(name), original_variable_(original_variable),
      device_shared_data_field_(nullptr){};
//=================================================================================================//
template <typename DataType>
void DiscreteDeviceSharedVariable<DataType>::allocateDataField(const size_t size)
{
    device_shared_data_field_ = allocateDeviceShared(size);
    original_variable_.setDelegateValue(device_shared_data_field_);
}
//=================================================================================================//
template <typename DataType>
DiscreteDeviceSharedVariable<DataType>::~DiscreteDeviceSharedVariable()
{
    freeDeviceData(device_shared_data_field_);
}
//=================================================================================================//
} // namespace SPH

#endif // BASE_VARIABLE_SYCL_HPP