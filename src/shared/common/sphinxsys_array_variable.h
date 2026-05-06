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
    ArrayData(DataType *data, size_t array_size)
        : data_(data), array_size_(array_size) {};

    size_t ArraySize() { return array_size_; };

    DataType *operator[](size_t particle_index)
    {
        return data_ + particle_index * array_size_;
    }

  protected
    DataType *data_;
    UnsignedInt array_size_;
};

template <typename DataType>
class ArrayVariable : protected DiscreteVariable<DataType>
{
  public:
    ArrayVariable(const std::string &array_name, StdVec<std::string> names, UnsignedInt variable_size)
        : DiscreteVariable<DataType>(array_name, names.size() * variable_size),
          names_(names) {};

    template <class ExecutionPolicy>
    ArrayData<DataType> DelegatedArrayData(const ExecutionPolicy &ex_policy)
    {
        return ArrayData(DiscreteVariable<DataType>::DelegatedData(ex_policy), names_.size());
    };

  protected:
    StdVec<std::string> names_;

    UnsignedInt getVariableIndex(const std::string &variable_name)
    {
        for (UnsignedInt i = 0; i < names_.size(); ++i)
        {
            if (names_[i] == variable_name)
            {
                return i;
            }
        }
        std::cout << "\n Error: variable name '" << variable_name
                  << "' is not found in ArrayVariable '" << this->Name() << "'!" << std::endl;
        exit(1);
    };
};
} // namespace SPH
#endif // SPHINXSYS_ARRAY_VARIABLE_H
