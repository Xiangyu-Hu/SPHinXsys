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

namespace SPH
{
template <typename DataType>
class ArrayData
{
  public:
    ArrayData(DataType *data_array, size_t array_size)
        : data_array_(data_array), array_size_(array_size) {};

    size_t getArraySize() { return array_size_; };

    DataType *operator[](size_t particle_index)
    {
        return data_array_ + particle_index * array_size_;
    }

  protected
    DataType *data_array_;
    UnsignedInt array_size_;
};

template <typename DataType>
class ArrayVariable : protected DiscreteVariable<DataType>
{
  public:
    ArrayVariable(StdVec<DiscreteVariable<DataType> *> variables, const std::string &collective_name)
        : DiscreteVariable<DataType>(collective_name, variables[0].getDataSize() * variables.size()),
          variables_(variables), array_size_(variables_.size()) {};

    template <class ExecutionPolicy>
    ArrayData<DataType> DelegatedArrayData(const ExecutionPolicy &ex_policy)
    {
        return ArrayData(DiscreteVariable<DataType>::DelegatedData(ex_policy), array_size_);
    };

  protected:
    StdVec<VariableType<DataType> *> variables_;
    UnsignedInt array_size_;
};
} // namespace SPH
#endif // SPHINXSYS_ARRAY_VARIABLE_H
