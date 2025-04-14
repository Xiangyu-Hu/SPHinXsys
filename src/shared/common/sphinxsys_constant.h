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
template <typename GeneratorType, typename DataType>
class ConstantArray : public Entity
{
    UniquePtrKeeper<Entity> device_only_constant_array_keeper_;

  public:
    ConstantArray(StdVec<GeneratorType *> generators)
        : Entity("ConstantArray"), generators_(generators),
          data_size_(generators.size()),
          data_(new DataType[data_size_]), delegated_(data_)
    {
        for (size_t i = 0; i != data_size_; ++i)
        {
            data_[i] = DataType(*generators_[i]);
        }
    };
    ~ConstantArray() { delete[] data_; };
    StdVec<GeneratorType *> getGenerators() { return generators_; }
    size_t getDataSize() { return data_size_; }
    DataType *Data() { return data_; };
    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return delegated_; };
    template <class PolicyType>
    DataType *DelegatedOnDevice(const DeviceExecution<PolicyType> &ex_policy);
    template <class PolicyType>
    DataType *DelegatedData(const DeviceExecution<PolicyType> &ex_policy) { return DelegatedOnDevice(ex_policy); };
    bool isDataDelegated() { return data_ != delegated_; };
    void setDelegateData(DataType *new_delegated) { delegated_ = new_delegated; };

  protected:
    StdVec<GeneratorType *> generators_;
    size_t data_size_;
    DataType *data_;
    DataType *delegated_;
};

template <typename GeneratorType, typename DataType>
class DeviceOnlyConstantArray : public Entity
{
  public:
    template <class PolicyType>
    DeviceOnlyConstantArray(const DeviceExecution<PolicyType> &ex_policy,
                            ConstantArray<GeneratorType, DataType> *host_constant);
    ~DeviceOnlyConstantArray();

  protected:
    DataType *device_only_data_;
};
} // namespace SPH
#endif // SPHINXSYS_CONSTANT_H
