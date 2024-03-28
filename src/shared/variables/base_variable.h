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
 * @file 	base_variable.h
 * @brief 	Here gives classes for the base variables used in simulation.
 * @details These variables are those discretized in spaces and time.
 * @author	Xiangyu Hu
 */

#ifndef BASE_VARIABLES_H
#define BASE_VARIABLES_H

#include "base_data_package.h"
namespace SPH
{
class BaseVariable
{
  public:
    explicit BaseVariable(const std::string &name) : name_(name){};
    virtual ~BaseVariable(){};
    std::string Name() const { return name_; };

  private:
    const std::string name_;
};

template <typename DataType>
using AllocatedMemory = DataType *;

template <typename DataType>
class Variable : public BaseVariable
{
    struct ConstantInitialization
    {
        explicit ConstantInitialization(DataType constant_value)
            : constant_value_(constant_value){};
        DataType operator()(size_t i) { return constant_value; };
    };

  public:
    explicit Variable(const std::string &name, size_t size = 1)
        : BaseVariable(name), size_(size),
          value_addrs_(std::malloc(size * std::sizeof(DataType))){};
    virtual ~Variable(){std::free(value_addrs_)};

    DataType *ValueAddress() { return value_addrs_; };

    template <typename InitializationFunction>
    void initializeVariable(const InitializationFunction &initialization)
    {
        for (size_t i = 0; i < size_; ++i)
        {
            value_addrs_[i] = initialization(i);
        }
    };

    void initializeVariable(DataType constant_value)
    {
        initializeVariable(ConstantInitialization(constant_value));
    };

  private:
    size_t size_;
    DataType *value_addrs_;
};

template <typename DataType>
Variable<DataType> *findVariableByName(DataContainerUniquePtrAssemble<Variable> &ptr_assemble, const std::string &name)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    auto &variables_ptrs = std::get<type_index>(ptr_assemble);
    auto result = std::find_if(variables.begin(), variables.end(),
                               [&](auto &variable) -> bool
                               { return variable.get()->Name() == name; });

    return result != variables.end() ? *result : nullptr;
};

template <typename DataType, typename... Args>
Variable<DataType> *addVariableToAssemble(DataContainerAssemble<AllocatedMemory> &assemble,
                                          DataContainerUniquePtrAssemble<Variable> &ptr_assemble, Args &&...args)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    auto &variable_ptrs = std::get<type_index>(ptr_assemble);
    Variable<DataType> *new_variable =
        variable_ptrs.template createPtr<VariableType<DataType>>(std::forward<Args>(args)...);
    std::get<type_index>(assemble).push_back(new_variable->ValueAddress());
    return new_variable;
};

template <typename DataType>
class DiscreteVariable : public BaseVariable
{
  public:
    DiscreteVariable(const std::string &name, std::size_t size)
        : BaseVariable(name), index_in_container_(index){};
    virtual ~DiscreteVariable(){};

    size_t IndexInContainer() const { return index_in_container_; };

  private:
    size_t index_in_container_;
};

template <typename DataType>
class MeshVariable : public BaseVariable
{
  public:
    MeshVariable(const std::string &name, size_t index)
        : BaseVariable(name), index_in_container_(index){};
    virtual ~MeshVariable(){};

    size_t IndexInContainer() const { return index_in_container_; };

  private:
    size_t index_in_container_;
};

template <typename DataType, template <typename VariableDataType> class VariableType>
VariableType<DataType> *findVariableByName(DataContainerAddressAssemble<VariableType> &assemble,
                                           const std::string &name)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    auto &variables = std::get<type_index>(assemble);
    auto result = std::find_if(variables.begin(), variables.end(),
                               [&](auto &variable) -> bool
                               { return variable->Name() == name; });

    return result != variables.end() ? *result : nullptr;
};

template <typename DataType, template <typename VariableDataType> class VariableType, typename... Args>
VariableType<DataType> *addVariableToAssemble(DataContainerAddressAssemble<VariableType> &assemble,
                                              DataContainerUniquePtrAssemble<VariableType> &ptr_assemble, Args &&...args)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    UniquePtrsKeeper<VariableType<DataType>> &variable_ptrs = std::get<type_index>(ptr_assemble);
    VariableType<DataType> *new_variable =
        variable_ptrs.template createPtr<VariableType<DataType>>(std::forward<Args>(args)...);
    std::get<type_index>(assemble).push_back(new_variable);
    return new_variable;
};
} // namespace SPH
#endif // BASE_VARIABLES_H
