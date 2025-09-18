#ifndef BASE_PARTICLES_HPP
#define BASE_PARTICLES_HPP

#include "base_particles.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::initializeVariable(DiscreteVariable<DataType> *variable, DataType initial_value)
{
    DataType *data_field = variable->Data();
    for (size_t i = 0; i != variable->getDataSize(); ++i)
    {
        data_field[i] = initial_value;
    }
    return data_field;
}
//=================================================================================================//
template <typename DataType, class InitializationFunction>
DataType *BaseParticles::
    initializeVariable(DiscreteVariable<DataType> *variable, const InitializationFunction &initialization)
{
    DataType *data_field = initializeVariable(variable);
    for (size_t i = 0; i != variable->getDataSize(); ++i)
    {
        data_field[i] = initialization(i); // Here, function object is applied for initialization.
    }
    return data_field;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::initializeVariable(
    DiscreteVariable<DataType> *variable, DiscreteVariable<DataType> *old_variable)
{
    DataType *data_field = variable->Data();
    DataType *old_data_field = old_variable->Data();
    for (size_t i = 0; i != variable->getDataSize(); ++i)
    {
        data_field[i] = old_data_field[i];
    }
    return data_field;
}
//=================================================================================================//
template <class DataType, typename... Args>
DataType *BaseParticles::
    addUniqueDiscreteVariableData(const std::string &name, size_t data_size, Args &&...args)
{
    DiscreteVariable<DataType> *variable =
        unique_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
    initializeVariable(variable, std::forward<Args>(args)...);
    return variable->Data();
}
//=================================================================================================//
template <class DataType, typename... Args>
DiscreteVariable<DataType> *BaseParticles::
    addUniqueDiscreteVariable(const std::string &name, size_t data_size, Args &&...args)
{
    DiscreteVariable<DataType> *variable =
        unique_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
    initializeVariable(variable, std::forward<Args>(args)...);
    return variable;
}
//=================================================================================================//
template <class DataType>
DiscreteVariable<DataType> *BaseParticles::addUniqueDiscreteVariableFrom(
    const std::string &name, DiscreteVariable<DataType> *old_variable)
{
    DataType *old_data_field = old_variable->Data();
    DiscreteVariable<DataType> *variable =
        unique_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, old_variable->getDataSize());
    initializeVariable(variable, [&](size_t index)
                       { return old_data_field[index]; });
    return variable;
}
//=================================================================================================//
template <typename DataType, typename... Args>
DataType *BaseParticles::registerDiscreteVariableData(const std::string &name,
                                                      size_t data_size, Args &&...args)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);
    if (variable == nullptr)
    {
        variable = addVariableToAssemble<DataType>(all_discrete_variables_, all_discrete_variable_ptrs_,
                                                   name, data_size);
        initializeVariable(variable, std::forward<Args>(args)...);
    }
    return variable->Data();
}
//=================================================================================================//
template <typename DataType, typename... Args>
DiscreteVariable<DataType> *BaseParticles::
    registerDiscreteVariable(const std::string &name, size_t data_size, Args &&...args)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);
    if (variable == nullptr)
    {
        variable = addVariableToAssemble<DataType>(all_discrete_variables_, all_discrete_variable_ptrs_,
                                                   name, data_size);
        initializeVariable(variable, std::forward<Args>(args)...);
    }
    return variable;
}
//=================================================================================================//
template <typename DataType, typename... Args>
DiscreteVariable<DataType> *BaseParticles::
    registerStateVariable(const std::string &name, Args &&...args)
{
    static_assert(DataTypeIndex<DataType>::value != DataTypeIndex<UnsignedInt>::value,
                  "\n Error: the data type UnsignedInt is not particle state variable!\n");

    DiscreteVariable<DataType> *variable =
        registerDiscreteVariable<DataType>(name, particles_bound_, std::forward<Args>(args)...);

    DataType *data_field = variable->Data();
    constexpr int type_index = DataTypeIndex<DataType>::value;
    std::get<type_index>(all_state_data_).push_back(data_field);

    return variable;
}
//=================================================================================================//
template <class DataType, typename... Args>
DataType *BaseParticles::addUniqueStateVariableData(const std::string &name, Args &&...args)
{
    return addUniqueDiscreteVariableData<DataType>(name, particles_bound_, std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename DataType, typename... Args>
DataType *BaseParticles::registerStateVariableData(const std::string &name, Args &&...args)
{

    constexpr int type_index = DataTypeIndex<DataType>::value;
    static_assert(DataTypeIndex<DataType>::value != DataTypeIndex<UnsignedInt>::value,
                  "\n Error: the data type UnsignedInt is not particle state variable!\n");

    DataType *data_field =
        registerDiscreteVariableData<DataType>(name, particles_bound_, std::forward<Args>(args)...);

    std::get<type_index>(all_state_data_).push_back(data_field);

    return data_field;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerStateVariableDataFrom(
    const std::string &new_name, const std::string &old_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, old_name);

    if (variable == nullptr)
    {
        std::cout << "\nError: the old variable '" << old_name << "' is not registered!\n";
        printBodyName();
        exit(1);
    }

    DataType *old_data_field = variable->Data();
    return registerStateVariableData<DataType>(new_name, [&](size_t index)
                                               { return old_data_field[index]; });
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerStateVariableDataFrom(
    const std::string &name, const StdVec<DataType> &geometric_data)
{
    DataType *data_field = registerStateVariableData<DataType>(name);

    for (size_t i = 0; i != geometric_data.size(); ++i)
    {
        data_field[i] = geometric_data[i];
    }
    return data_field;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerStateVariableDataFromReload(const std::string &name)
{
    DataType *data_field = registerStateVariableData<DataType>(name);

    size_t index = 0;
    for (auto child = reload_xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
    {
        reload_xml_parser_.queryAttributeValue(child, name, data_field[index]);
        index++;
    }

    return data_field;
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::registerStateVariableFrom(
    const std::string &new_name, const std::string &old_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, old_name);

    if (variable == nullptr)
    {
        std::cout << "\nError: the old variable '" << old_name << "' is not registered!\n";
        printBodyName();
        exit(1);
    }

    DataType *old_data_field = variable->Data();
    return registerStateVariable<DataType>(new_name, [&](size_t index)
                                           { return old_data_field[index]; });
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::registerStateVariableFrom(
    const std::string &name, const StdVec<DataType> &geometric_data)
{
    DiscreteVariable<DataType> *variable = registerStateVariable<DataType>(name);
    DataType *data_field = variable->Data();
    for (size_t i = 0; i != geometric_data.size(); ++i)
    {
        data_field[i] = geometric_data[i];
    }
    return variable;
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::registerStateVariableFromReload(const std::string &name)
{
    DiscreteVariable<DataType> *new_variable = registerStateVariable<DataType>(name);
    DataType *data_field = new_variable->Data();

    size_t index = 0;
    for (auto child = reload_xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
    {
        reload_xml_parser_.queryAttributeValue(child, name, data_field[index]);
        index++;
    }

    return new_variable;
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::getVariableByName(const std::string &name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);
    if (variable == nullptr)
    {
        std::cout << "\nError: the variable '" << name << "' is not registered!\n";
        printBodyName();
        exit(1);
    }
    return variable;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getVariableDataByName(const std::string &name)
{
    DiscreteVariable<DataType> *variable = getVariableByName<DataType>(name);

    if (variable->Data() == nullptr)
    {
        std::cout << "\nError: the variable '" << name << "' has not been allocated yet!\n";
        printBodyName();
        exit(1);
    }

    return variable->Data();
}
//=================================================================================================//
template <typename DataType>
StdVec<DiscreteVariable<DataType> *> BaseParticles::registerStateVariables(
    const StdVec<std::string> &names, const std::string &suffix)
{
    StdVec<DiscreteVariable<DataType> *> variables;
    for (auto &name : names)
    {
        std::string variable_name = name + suffix;
        variables.push_back(registerStateVariable<DataType>(variable_name));
    }
    return variables;
}
//=================================================================================================//
template <typename DataType>
StdVec<DiscreteVariable<DataType> *> BaseParticles::getVariablesByName(
    const StdVec<std::string> &names, const std::string &suffix)
{
    StdVec<DiscreteVariable<DataType> *> variables;
    for (auto &name : names)
    {
        std::string variable_name = name + suffix;
        variables.push_back(getVariableByName<DataType>(variable_name));
    }
    return variables;
}
//=================================================================================================//
template <class DataType>
SingularVariable<DataType> *BaseParticles::
    addUniqueSingularVariable(const std::string &name, DataType initial_value)
{

    SingularVariable<DataType> *variable =
        unique_variable_ptrs_.createPtr<SingularVariable<DataType>>(name, initial_value);
    return variable;
}
//=================================================================================================//
template <typename DataType>
SingularVariable<DataType> *BaseParticles::
    registerSingularVariable(const std::string &name, DataType initial_value)
{
    SingularVariable<DataType> *variable = findVariableByName<DataType>(all_singular_variables_, name);

    return variable != nullptr
               ? variable
               : addVariableToAssemble<DataType>(
                     all_singular_variables_, all_singular_variable_ptrs_, name, initial_value);
}
//=================================================================================================//
template <typename DataType>
SingularVariable<DataType> *BaseParticles::getSingularVariableByName(const std::string &name)
{
    SingularVariable<DataType> *variable = findVariableByName<DataType>(all_singular_variables_, name);

    if (variable == nullptr)
    {
        std::cout << "\nError: the variable '" << name << "' is not registered!\n";
        printBodyName();
    }

    return variable;
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::
    addVariableToList(ParticleVariables &variable_set, const std::string &name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);

    if (variable == nullptr)
    {
        std::cout << "\n Error: the " << type_name<DataType>() << " variable '" << name << "' is  not exist!" << std::endl;
        printBodyName();
        exit(1);
    }

    return addVariableToList<DataType>(variable_set, variable);
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::
    addVariableToList(ParticleVariables &variable_set, DiscreteVariable<DataType> *variable)
{
    if (variable->getDataSize() < particles_bound_)
    {
        std::cout << "\n Error: The variable '" << variable->Name() << "' can not be treated as a particle variable," << std::endl;
        std::cout << "\n because the data size " << variable->getDataSize() << " is too less!" << std::endl;
        printBodyName();
        exit(1);
    }

    DiscreteVariable<DataType> *listed_variable = findVariableByName<DataType>(variable_set, variable->Name());
    if (listed_variable == nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        std::get<type_index>(variable_set).push_back(variable);
        return variable;
    }

    return nullptr; // no variable added as sortable variable
}
//=================================================================================================//
template <typename DataType, typename... Args>
void BaseParticles::addEvolvingVariable(Args &&...args)
{
    DiscreteVariable<DataType> *new_sortable =
        addVariableToList<DataType>(evolving_variables_, std::forward<Args>(args)...);
    if (new_sortable != nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        DataType *data_field = new_sortable->Data();
        std::get<type_index>(evolving_variables_data_).push_back(data_field);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addEvolvingVariable(DiscreteVariableArray<DataType> *variable_array)
{
    StdVec<DiscreteVariable<DataType> *> variables = variable_array->getVariables();
    for (size_t i = 0; i != variables.size(); ++i)
    {
        addEvolvingVariable<DataType>(variables[i]);
    }
}
//=================================================================================================//
template <typename DataType, typename... Args>
void BaseParticles::addVariableToWrite(Args &&...args)
{
    addVariableToList<DataType>(variables_to_write_, std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToWrite(DiscreteVariableArray<DataType> *variable_array)
{
    StdVec<DiscreteVariable<DataType> *> variables = variable_array->getVariables();
    for (size_t i = 0; i != variables.size(); ++i)
    {
        addVariableToWrite<DataType>(variables[i]);
    }
}
//===============================================================================
template <typename DataType>
BaseParticles *BaseParticles::reloadExtraVariable(const std::string &name)
{
    registerStateVariableFromReload<DataType>(name);
    return this;
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::CopyParticleState::
operator()(DataContainerKeeper<AllocatedData<DataType>> &data_keeper, size_t index, size_t another_index)
{
    for (size_t i = 0; i != data_keeper.size(); ++i)
    {
        data_keeper[i][index] = data_keeper[i][another_index];
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::WriteAParticleVariableToXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, XmlParser &xml_parser)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        DataType *data_field = variables[i]->Data();
        for (auto child = xml_parser.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser.setAttributeToElement(child, variables[i]->Name(), data_field[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ReadAParticleVariableFromXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
           BaseParticles *base_particles, XmlParser &xml_parser)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        DataType *data_field = variables[i]->Data() != nullptr
                                   ? variables[i]->Data()
                                   : base_particles->initializeVariable<DataType>(variables[i]);
        for (auto child = xml_parser.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser.queryAttributeValue(child, variables[i]->Name(), data_field[index]);
            index++;
        }
    }
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLES_HPP
