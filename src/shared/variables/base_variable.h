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
#include <cstring>
#include <stdio.h>
#include <typeinfo>

namespace SPH
{
class Entity
{
  public:
    explicit Entity(const std::string &name)
        : name_(name), entity_id_(0){};
    virtual ~Entity(){};
    std::string Name() const { return name_; };
    size_t EntityID() const
    {
        if (entity_id_ == 0)
        {
            std::cout << "Entity ID is not set for " << name_ << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        return entity_id_;
    };
    void setEntityID(size_t entity_id) { entity_id_ = entity_id; };

  private:
    const std::string name_; // name of an entity object
    size_t entity_id_;       // unique id for each derived entity type
};

template <typename FunctionalType>
class FunctionalEntity : public Entity
{
  public:
    template <typename... Args>
    explicit FunctionalEntity(const std::string &name, Args &&...args)
        : Entity(name), functional_(new FunctionalType(std::forward<Args>(args)...))
    {
        setEntityID(typeid(FunctionalType).hash_code());
    };
    virtual ~FunctionalEntity() { delete functional_; };
    FunctionalType *Functional() { return functional_; };

  private:
    FunctionalType *functional_;
};

template <typename DataType>
class SingularVariable : public Entity
{
  public:
    SingularVariable(const std::string &name, DataType &value)
        : Entity(name), value_(value){};
    virtual ~SingularVariable(){};

    DataType *ValueAddress() { return &value_; };

  private:
    DataType value_;
};

template <typename DataType>
class DiscreteVariable : public Entity
{
  public:
    DiscreteVariable(const std::string &name)
        : Entity(name), data_field_(nullptr){};
    virtual ~DiscreteVariable() { delete[] data_field_; };
    DataType *DataField() { return data_field_; };
    void allocateDataField(const size_t size)
    {
        data_field_ = new DataType[size];
    }

  private:
    DataType *data_field_;
};

template <typename DataType>
class MeshVariable : public Entity
{
  public:
    using PackageData = PackageDataMatrix<DataType, 4>;
    MeshVariable(const std::string &name, size_t data_size)
        : Entity(name), data_field_(nullptr){};
    virtual ~MeshVariable() { delete[] data_field_; };

    // void setDataField(PackageData* mesh_data){ data_field_ = mesh_data; };
    PackageData *DataField() { return data_field_; };
    void allocateAllMeshVariableData(const size_t size)
    {
        data_field_ = new PackageData[size];
    }

  private:
    PackageData *data_field_;
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
