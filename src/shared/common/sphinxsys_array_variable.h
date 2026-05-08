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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file sphinxsys_array_variable.h
 * @brief tbd
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_ARRAY_VARIABLE_H
#define SPHINXSYS_ARRAY_VARIABLE_H

#include "sphinxsys_variable.h"
#include "sphinxsys_variable_array.h"
namespace SPH
{
template <typename DataType>
class ArrayData // transposed view of DataArray
{
  public:
    ArrayData() : transposed_data_(nullptr), array_size_(0) {};
    ArrayData(DataType *transposed_data, size_t array_size)
        : transposed_data_(transposed_data), array_size_(array_size) {};

    size_t ArraySize() { return array_size_; };

    DataType *operator[](size_t particle_index)
    {
        return transposed_data_ + particle_index * array_size_;
    }

  protected:
    DataType *transposed_data_;
    UnsignedInt array_size_;
};

template <typename DataType>
class ArrayVariable : protected DiscreteVariable<DataType>
{
  public:
    ArrayVariable(const std::string &name, UnsignedInt variable_size, UnsignedInt array_size)
        : DiscreteVariable<DataType>(name, variable_size * array_size);

    template <class ExecutionPolicy>
    ArrayData<DataType> DelegatedArrayData(const ExecutionPolicy &ex_policy)
    {
        return ArrayData(DiscreteVariable<DataType>::DelegatedData(ex_policy), names_.size());
    };
};

/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<ArrayVariable> ArrayVariables;
} // namespace SPH
#endif // SPHINXSYS_ARRAY_VARIABLE_H
