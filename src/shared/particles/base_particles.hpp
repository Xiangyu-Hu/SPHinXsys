/**
 * @file 	base_particles.hpp
 * @brief 	This is the implementation of the template functions in base_particles.h
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_PARTICLES_HPP
#define BASE_PARTICLES_HPP

#include "base_particles.h"
#include "particle_dynamics_algorithms.h"

//=====================================================================================================//
namespace SPH
{
//=================================================================================================//
template <typename DataType>
void BaseParticles::registerVariable(StdLargeVec<DataType> &variable_addrs,
                                     const std::string &variable_name, DataType initial_value)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    if (variable == nullptr)
    {
        variable_addrs.resize(particles_bound_, initial_value);

        constexpr int type_index = DataTypeIndex<DataType>::value;
        std::get<type_index>(all_particle_data_).push_back(&variable_addrs);
        size_t new_variable_index = std::get<type_index>(all_particle_data_).size() - 1;

        addVariableToAssemble<DataType>(all_discrete_variables_, all_discrete_variable_ptrs_, variable_name, new_variable_index);
    }
    else
    {
        std::cout << "\n Error: the variable '" << variable_name << "' has already been registered!" << std::endl;
        std::cout << "\n Please check if " << variable_name << " is a sharable variable." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
template <typename DataType, class InitializationFunction>
void BaseParticles::registerVariable(StdLargeVec<DataType> &variable_addrs,
                                     const std::string &variable_name, const InitializationFunction &initialization)
{
    registerVariable(variable_addrs, variable_name);
    for (size_t i = 0; i != particles_bound_; ++i)
    {
        variable_addrs[i] = initialization(i); // Here, lambda function is applied for initialization.
    }
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerSingleVariable(const std::string &variable_name, DataType initial_value)
{
    SingleVariable<DataType> *variable = findVariableByName<DataType>(all_single_variables_, variable_name);

    return variable != nullptr
               ? variable->ValueAddress()
               : addVariableToAssemble<DataType>(all_single_variables_,
                                                 all_global_variable_ptrs_, variable_name, initial_value)
                     ->ValueAddress();
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getSingleVariableByName(const std::string &variable_name)
{
    SingleVariable<DataType> *variable = findVariableByName<DataType>(all_single_variables_, variable_name);

    if (variable != nullptr)
    {
        return variable->ValueAddress();
    }

    std::cout << "\nError: the variable '" << variable_name << "' is not registered!\n";
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    return nullptr;

    return nullptr;
}
//=================================================================================================//
template <typename DataType, typename... Args>
StdLargeVec<DataType> *BaseParticles::registerSharedVariable(const std::string &variable_name, Args &&...args)
{

    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    constexpr int type_index = DataTypeIndex<DataType>::value;
    if (variable == nullptr)
    {
        UniquePtrsKeeper<StdLargeVec<DataType>> &container = std::get<type_index>(shared_particle_data_ptrs_);
        StdLargeVec<DataType> *contained_data = container.template createPtr<StdLargeVec<DataType>>();
        registerVariable(*contained_data, variable_name, std::forward<Args>(args)...);
        return contained_data;
    }
    else
    {
        return std::get<type_index>(all_particle_data_)[variable->IndexInContainer()];
    }
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::registerSharedVariableFrom(const std::string &new_name, const std::string &old_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, old_name);

    if (variable != nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        StdLargeVec<DataType> *old_data = std::get<type_index>(all_particle_data_)[variable->IndexInContainer()];
        return registerSharedVariable<DataType>(new_name, [&](size_t index)
                                                { return (*old_data)[index]; });
    }
    else
    {
        std::cout << "\nError: the old variable '" << old_name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    return nullptr;
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::getVariableByName(const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    if (variable != nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        return std::get<type_index>(all_particle_data_)[variable->IndexInContainer()];
    }
    else
    {
        std::cout << "\nError: the variable '" << variable_name << "' is not registered!\n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    return nullptr;
}
//=================================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *BaseParticles::
    addVariableToList(ParticleVariables &variable_set, const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    if (variable != nullptr)
    {
        DiscreteVariable<DataType> *listed_variable = findVariableByName<DataType>(variable_set, variable_name);

        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(variable_set).push_back(variable);
            return variable;
        }
    }
    else
    {
        std::cout << "\n Error: the variable '" << variable_name << "' to write is not particle data!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return nullptr;
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToSort(const std::string &variable_name)
{
    DiscreteVariable<DataType> *new_sortable =
        addVariableToList<DataType>(sortable_variables_, variable_name);
    if (new_sortable != nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        StdLargeVec<DataType> *variable_data = std::get<type_index>(all_particle_data_)[new_sortable->IndexInContainer()];
        std::get<type_index>(sortable_data_).push_back(variable_data);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToWrite(const std::string &variable_name)
{
    addVariableToList<DataType>(variables_to_write_, variable_name);
}
//=================================================================================================//
template <class DerivedVariableMethod, class... Ts>
void BaseParticles::addDerivedVariableToWrite(SPHBody &sph_body, Ts &&...args)
{
    SimpleDynamics<DerivedVariableMethod> *derived_data =
        derived_particle_data_.createPtr<SimpleDynamics<DerivedVariableMethod>>(sph_body, std::forward<Ts>(args)...);
    derived_variables_.push_back(derived_data);
    using DerivedDataType = typename DerivedVariableMethod::DerivedDataType;
    addVariableToList<DerivedDataType>(variables_to_write_, derived_data->variable_name_);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToRestart(const std::string &variable_name)
{
    addVariableToList<DataType>(variables_to_restart_, variable_name);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToReload(const std::string &variable_name)
{
    addVariableToList<DataType>(variables_to_reload_, variable_name);
}
//=================================================================================================//
template <typename SequenceMethod>
void BaseParticles::sortParticles(SequenceMethod &sequence_method)
{
    StdLargeVec<size_t> &sequence = sequence_method.computingSequence(*this);
    particle_sorting_.sortingParticleData(sequence.data(), total_real_particles_);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ResizeParticles::
operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper, size_t new_size)
{
    for (size_t i = 0; i != data_keeper.size(); ++i)
    {
        data_keeper[i]->resize(new_size, ZeroData<DataType>::value);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::CopyParticleData::
operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper, size_t index, size_t another_index)
{
    for (size_t i = 0; i != data_keeper.size(); ++i)
    {
        (*data_keeper[i])[index] = (*data_keeper[i])[another_index];
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::WriteAParticleVariableToXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, ParticleData &all_particle_data)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        StdLargeVec<DataType> &variable_data = *(std::get<type_index>(all_particle_data)[variables[i]->IndexInContainer()]);
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.setAttributeToElement(child, variables[i]->Name(), variable_data[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ReadAParticleVariableFromXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, ParticleData &all_particle_data)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        StdLargeVec<DataType> &variable_data = *(std::get<type_index>(all_particle_data)[variables[i]->IndexInContainer()]);
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.queryAttributeValue(child, variables[i]->Name(), variable_data[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename StreamType>
void BaseParticles::writeParticlesToVtk(StreamType &output_stream)
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

    // write unsorted particles ID
    output_stream << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        output_stream << unsorted_id_[i] << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";

    // write integers
    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        StdLargeVec<int> &variable_data = *(std::get<type_index_int>(all_particle_data_)[variable->IndexInContainer()]);
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Int32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << variable_data[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write scalars
    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        StdLargeVec<Real> &variable_data = *(std::get<type_index_Real>(all_particle_data_)[variable->IndexInContainer()]);
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << variable_data[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write vectors
    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        StdLargeVec<Vecd> &variable_data = *(std::get<type_index_Vecd>(all_particle_data_)[variable->IndexInContainer()]);
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Vec3d vector_value = upgradeToVec3d(variable_data[i]);
            output_stream << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write matrices
    constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
    for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(variables_to_write_))
    {
        StdLargeVec<Matd> &variable_data = *(std::get<type_index_Matd>(all_particle_data_)[variable->IndexInContainer()]);
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type= \"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Mat3d matrix_value = upgradeToMat3d(variable_data[i]);
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
template <typename DataType>
BaseDerivedVariable<DataType>::
    BaseDerivedVariable(SPHBody &sph_body, const std::string &variable_name)
    : variable_name_(variable_name)
{
    sph_body.getBaseParticles().registerVariable(derived_variable_, variable_name_);
};
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLES_HPP