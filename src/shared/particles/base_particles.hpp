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
template <typename DataType>
DataType *BaseParticles::allocateVariable(DiscreteVariable<DataType> *variable)
{
    if (variable->DataField() == nullptr)
    {
        variable->allocateDataField(particles_bound_);
    }
    else
    {
        std::cout << "\n Error: the variable '" << variable->Name() << "' has already been allocated!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return variable->DataField();
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::initializeVariable(DiscreteVariable<DataType> *variable, DataType initial_value)
{
    DataType *data_field = allocateVariable(variable);
    for (size_t i = 0; i != particles_bound_; ++i)
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
    for (size_t i = 0; i != total_real_particles_; ++i)
    {
        data_field[i] = initialization(i); // Here, function object is applied for initialization.
    }
    return data_field;
}
//=================================================================================================//
template <typename DataType, typename... Args>
DataType *BaseParticles::addUniqueDiscreteVariable(const std::string &name, Args &&...args)
{

    DiscreteVariable<DataType> *variable =
        unique_discrete_variable_ptrs_.createPtr<DiscreteVariable<DataType>>(name);
    initializeVariable(variable, std::forward<Args>(args)...);
    return variable->DataField();
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerSingularVariable(const std::string &name, DataType initial_value)
{
    SingularVariable<DataType> *variable = findVariableByName<DataType>(all_singular_variables_, name);

    return variable != nullptr
               ? variable->ValueAddress()
               : addVariableToAssemble<DataType>(all_singular_variables_,
                                                 all_global_variable_ptrs_, name, initial_value)
                     ->ValueAddress();
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
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::addSharedVariable(const std::string &name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);
    if (variable == nullptr)
    {
        variable = addVariableToAssemble<DataType>(all_discrete_variables_, all_discrete_variable_ptrs_, name);
    }
    return variable;
}
//=================================================================================================//
template <typename DataType, typename... Args>
DataType *BaseParticles::registerSharedVariable(const std::string &name, Args &&...args)
{

    DiscreteVariable<DataType> *variable = addSharedVariable<DataType>(name);

    if (variable->DataField() == nullptr)
    {
        initializeVariable(variable, std::forward<Args>(args)...);
        constexpr int type_index = DataTypeIndex<DataType>::value;
        if (type_index != DataTypeIndex<size_t>::value) // particle IDs excluded
        {
            std::get<type_index>(all_state_data_).push_back(variable->DataField());
        }
    }

    return variable->DataField();
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerSharedVariableFrom(
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
    return registerSharedVariable<DataType>(new_name, [&](size_t index)
                                            { return old_data_field[index]; });
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerSharedVariableFrom(
    const std::string &name, const StdLargeVec<DataType> &geometric_data)
{
    return registerSharedVariable<DataType>(name, [&](size_t index)
                                            { return geometric_data[index]; });
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerSharedVariableFromReload(const std::string &name)
{
    DataType *data_field = registerSharedVariable<DataType>(name);

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
template <typename OutStreamType>
void BaseParticles::writeParticlesToVtk(OutStreamType &output_stream)
{
    size_t total_real_particles = total_real_particles_;

    // write sorted particles ID
    output_stream << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        output_stream << i << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";

    // write original particles ID
    output_stream << "    <DataArray Name=\"OriginalParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        output_stream << original_id_[i] << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";

    // write integers
    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        int *data_field = variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Int32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << data_field[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write scalars
    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        Real *data_field = variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << data_field[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write vectors
    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        Vecd *data_field = variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Vec3d vector_value = upgradeToVec3d(data_field[i]);
            output_stream << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write matrices
    constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
    for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(variables_to_write_))
    {
        Matd *data_field = variable->DataField();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type= \"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Mat3d matrix_value = upgradeToMat3d(data_field[i]);
            for (int k = 0; k != 3; ++k)
            {
                Vec3d col_vector = matrix_value.col(k);
                output_stream << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
            }
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLES_HPP