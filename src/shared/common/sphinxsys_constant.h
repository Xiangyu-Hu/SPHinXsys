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
template <typename GeneratorType, typename FunctionType>
class Constant : public Entity
{
    UniquePtrKeeper<Entity> device_shared_constant_keeper_;

  public:
    using DataType = typename GeneratorType::FunctionType;
    Constant(GeneratorType *generator, size_t data_size = 1)
        : Entity(generator->getFunctionName<DataType>()),
          generator_(generator), data_size_(data_size),
          data_(new DataType[data_size]), delegated_(data_)
    {
        for (size_t i = 0; i != data_size_; ++i)
        {
            data_[i] = generator_[i]->getFunction<DataType>();
        }
    };
    ~Constant() { delete[] data_; };
    GeneratorType *getGenerator() { return generator_; }
    size_t getDataSize() { return data_size_; }
    DataType *Data() { return data_; };
    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return delegated_; };
    DataType *DelegatedData(const ParallelDevicePolicy &par_device);
    bool isDataDelegated() { return data_ != delegated_; };
    void setDelegateData(DataType *new_delegated) { delegated_ = new_delegated; };

  protected:
    GeneratorType *generator_;
    size_t data_size_;
    DataType *data_;
    DataType *delegated_;
};

template <typename GeneratorType, typename FunctionType>
class DeviceOnlyConstant : public Entity
{
    using DataType = typename GeneratorType::FunctionType;

  public:
    DeviceOnlyConstant(Constant<GeneratorType, FunctionType> *host_constant);
    ~DeviceOnlyConstant();

  protected:
    DataType *device_only_data_;
};
} // namespace SPH
#endif // SPHINXSYS_CONSTANT_H
