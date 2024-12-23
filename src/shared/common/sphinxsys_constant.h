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
 * @file sphinxsys_constant.h
 * @brief Here gives classes for the constants and variables used in simulation.
 * @details These variables are those discretized in spaces and time.
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_CONSTANT_H
#define SPHINXSYS_CONSTANT_H

#include "sphinxsys_variable.h"

namespace SPH
{
template <typename DataType>
class DiscreteConstant;

template <typename DataType>
using SingularConstant = SingularVariable<DataType>;

template <typename DataType>
using DeviceSharedSingularConstant = DeviceSharedSingularVariable<DataType>;

template <typename DataType>
class DeviceSharedDiscreteConstant : public Entity
{
  public:
    DeviceSharedDiscreteConstant(DiscreteConstant<DataType> *host_variable);
    ~DeviceSharedDiscreteConstant();

  protected:
    DataType *device_shared_data_field_;
};

template <typename DataType>
class DiscreteConstant : public Entity
{
    UniquePtrKeeper<Entity> device_shared_constant_keeper_;

  public:
    DiscreteConstant(const std::string &name, size_t data_size)
        : Entity(name), data_size_(data_size), data_field_(new DataType[data_size]),
          delegated_data_field_(data_field_){};
    ~DiscreteConstant() { delete[] data_field_; };
    bool isDataDelegated() { return data_field_ != delegated_data_field_; };
    size_t getDataSize() { return data_size_; }
    DataType *Data() { return delegated_data_field_; };
    void setDeviceData(DataType *data_field) { delegated_data_field_ = data_field; };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return delegated_data_field_; };
    DataType *DelegatedData(const ParallelDevicePolicy &par_device)
    {
        if (!isDataDelegated())
        {
            device_shared_constant_keeper_
                .createPtr<DeviceSharedDiscreteConstant<DataType>>(this);
        }
        return delegated_data_field_;
    };

  private:
    size_t data_size_;
    DataType *data_field_;
    DataType *delegated_data_field_;
};
} // namespace SPH
#endif // SPHINXSYS_CONSTANT_H
