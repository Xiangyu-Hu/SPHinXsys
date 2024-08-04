#ifndef BASE_PARTICLES_HPP
#define BASE_PARTICLES_HPP

#include "base_particles.h"

namespace SPH
{
//=================================================================================================//
template <typename OwnerType>
void BaseParticles::checkReloadFileRead(OwnerType *owner)
{
    if (reload_xml_parser_.first_element_ == nullptr)
    {
        std::cout << "\n Error: the reload file is not read! \n";
        std::cout << "\n This error occurs in " << typeid(*owner).name() << '\n';
        exit(1);
    }
}
//=================================================================================================//
template <typename DataType, template <typename T> class VariableType>
DataType *BaseParticles::initializeVariable(VariableType<DataType> *variable, DataType initial_value)
{
    DataType *data_field = variable->DataField();
    for (size_t i = 0; i != variable->getDataSize(); ++i)
    {
        data_field[i] = initial_value;
    }
    return data_field;
}
//=================================================================================================//
template <typename DataType, template <typename T> class VariableType, class InitializationFunction>
DataType *BaseParticles::
    initializeVariable(VariableType<DataType> *variable, const InitializationFunction &initialization)
{
    DataType *data_field = initializeVariable(variable);
    for (size_t i = 0; i != variable->getDataSize(); ++i)
    {
        data_field[i] = initialization(i); // Here, function object is applied for initialization.
    }
    return data_field;
}
//=================================================================================================//
template <class DataType, typename... Args>
DataType *BaseParticles::
    addUniqueDiscreteVariable(const std::string &name, size_t data_size, Args &&...args)
{

    DiscreteVariable<DataType> *variable =
        unique_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name, data_size);
    initializeVariable(variable, std::forward<Args>(args)...);
    return variable->DataField();
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
                     all_singular_variables_, all_global_variable_ptrs_, name, initial_value);
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getSingularVariableByName(const std::string &name)
{
    SingularVariable<DataType> *variable = findVariableByName<DataType>(all_singular_variables_, name);

    if (variable != nullptr)
    {
        return variable->ValueAddress();
    }

    std::cout << "\nError: the variable '" << name << "' is not registered!\n";
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    return nullptr;

    return nullptr;
}
//=================================================================================================//
template <typename DataType, typename... Args>
DataType *BaseParticles::registerDiscreteVariable(const std::string &name,
                                                  size_t data_size, Args &&...args)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);
    if (variable == nullptr)
    {
        variable = addVariableToAssemble<DataType>(all_discrete_variables_, all_discrete_variable_ptrs_,
                                                   name, data_size);
        initializeVariable(variable, std::forward<Args>(args)...);
    }
    return variable->DataField();
}
//=================================================================================================//
template <typename DataType, typename... Args>
DataType *BaseParticles::registerStateVariable(const std::string &name, Args &&...args)
{

    constexpr int type_index = DataTypeIndex<DataType>::value;
    static_assert(DataTypeIndex<DataType>::value != DataTypeIndex<UnsignedInt>::value,
                  "\n Error: the data type UnsignedInt is not particle state variable!\n");

    DataType *data_field =
        registerDiscreteVariable<DataType>(name, particles_bound_, std::forward<Args>(args)...);

    std::get<type_index>(all_state_data_).push_back(data_field);

    return data_field;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerStateVariableFrom(
    const std::string &new_name, const std::string &old_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, old_name);

    if (variable == nullptr)
    {
        std::cout << "\nError: the old variable '" << old_name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    DataType *old_data_field = variable->DataField();
    return registerStateVariable<DataType>(new_name, [&](size_t index)
                                           { return old_data_field[index]; });
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerStateVariableFrom(
    const std::string &name, const StdLargeVec<DataType> &geometric_data)
{
    DataType *data_field = registerStateVariable<DataType>(name);

    for (size_t i = 0; i != geometric_data.size(); ++i)
    {
        data_field[i] = geometric_data[i];
    }
    return data_field;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerStateVariableFromReload(const std::string &name)
{
    DataType *data_field = registerStateVariable<DataType>(name);

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
DiscreteVariable<DataType> *BaseParticles::getVariableByName(const std::string &name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);
    if (variable == nullptr)
    {
        std::cout << "\nError: the variable '" << name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return variable;
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getVariableDataByName(const std::string &name)
{
    DiscreteVariable<DataType> *variable = getVariableByName<DataType>(name);

    if (variable->DataField() == nullptr)
    {
        std::cout << "\nError: the variable '" << name << "' has not been allocated yet!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    return variable->DataField();
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::
    addVariableToList(ParticleVariables &variable_set, const std::string &name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);

    if (variable != nullptr)
    {
        DiscreteVariable<DataType> *listed_variable = findVariableByName<DataType>(variable_set, name);

        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(variable_set).push_back(variable);
            return variable;
        }
    }
    else
    {
        std::cout << "\n Error: the variable '" << name << "' to write is not particle data!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return nullptr;
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToSort(const std::string &name)
{
    DiscreteVariable<DataType> *new_sortable =
        addVariableToList<DataType>(sortable_variables_, name);
    if (new_sortable != nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        DataType *data_field = new_sortable->DataField();
        std::get<type_index>(sortable_data_).push_back(data_field);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToWrite(const std::string &name)
{
    addVariableToList<DataType>(variables_to_write_, name);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToRestart(const std::string &name)
{
    addVariableToList<DataType>(variables_to_restart_, name);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToReload(const std::string &name)
{
    addVariableToList<DataType>(variables_to_reload_, name);
}
//=================================================================================================//
template <typename DataType, typename... Args>
DataType *BaseParticles::
    registerStateVariable(const ParallelDevicePolicy &execution_policy, const std::string &name, Args &&...args)
{
    registerStateVariable<DataType>(name, std::forward<Args>(args)...);
    return getVariableDataByName<DataType>(execution_policy, name);
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getVariableDataByName(const ParallelDevicePolicy &execution_policy, const std::string &name)
{
    DiscreteVariable<DataType> *variable = getVariableByName<DataType>(name);
    if (!variable->existDeviceDataField())
    {
        unique_variable_ptrs_.createPtr<DiscreteDeviceOnlyVariable<DataType>>(variable);
    }
    return variable->DeviceDataField();
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
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        DataType *data_field = variables[i]->DataField();
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.setAttributeToElement(child, variables[i]->Name(), data_field[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ReadAParticleVariableFromXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, BaseParticles *base_particles)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        DataType *data_field = variables[i]->DataField() != nullptr
                                   ? variables[i]->DataField()
                                   : base_particles->initializeVariable<DataType>(variables[i]);
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.queryAttributeValue(child, variables[i]->Name(), data_field[index]);
            index++;
        }
    }
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLES_HPP