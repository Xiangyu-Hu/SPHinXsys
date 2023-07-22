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
        variable_addrs.resize(real_particles_bound_, initial_value);

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
    for (size_t i = 0; i != real_particles_bound_; ++i)
    {
        variable_addrs[i] = initialization(i); // Here, lambda function is applied for initialization.
    }
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerGlobalVariable(const std::string &variable_name, DataType initial_value)
{
    GlobalVariable<DataType> *variable = findVariableByName<DataType>(all_global_variables_, variable_name);

    return variable != nullptr
               ? variable->ValueAddress()
               : addVariableToAssemble<DataType>(all_global_variables_,
                                                 all_global_variable_ptrs_, variable_name, initial_value)
                     ->ValueAddress();
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getGlobalVariableByName(const std::string &variable_name)
{
    GlobalVariable<DataType> *variable = findVariableByName<DataType>(all_global_variables_, variable_name);

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
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::
    registerSharedVariable(const std::string &variable_name, const DataType &default_value)
{

    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    constexpr int type_index = DataTypeIndex<DataType>::value;
    if (variable == nullptr)
    {
        UniquePtrsKeeper<StdLargeVec<DataType>> &container = std::get<type_index>(shared_particle_data_ptrs_);
        StdLargeVec<DataType> *contained_data = container.template createPtr<StdLargeVec<DataType>>();
        registerVariable(*contained_data, variable_name, default_value);
        return contained_data;
    }
    else
    {
        return std::get<type_index>(all_particle_data_)[variable->IndexInContainer()];
    }
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

    std::cout << "\nError: the variable '" << variable_name << "' is not registered!\n";
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    return nullptr;

    return nullptr;
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToList(ParticleVariables &variable_set, const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    if (variable != nullptr)
    {
        DiscreteVariable<DataType> *listed_variable = findVariableByName<DataType>(variable_set, variable_name);

        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(variable_set).push_back(variable);
        }
    }
    else
    {
        std::cout << "\n Error: the variable '" << variable_name << "' to write is not particle data!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
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
void BaseParticles::addDerivedVariableToWrite(Ts &&...args)
{
    SimpleDynamics<DerivedVariableMethod> *derived_data =
        derived_particle_data_.createPtr<SimpleDynamics<DerivedVariableMethod>>(sph_body_, std::forward<Ts>(args)...);
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
template <typename DataType>
void BaseParticles::registerSortableVariable(const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, variable_name);

    if (variable != nullptr)
    {
        DiscreteVariable<DataType> *listed_variable = findVariableByName<DataType>(sortable_variables_, variable_name);

        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(sortable_variables_).push_back(variable);
            StdLargeVec<DataType> *variable_data = std::get<type_index>(all_particle_data_)[variable->IndexInContainer()];
            std::get<type_index>(sortable_data_).push_back(variable_data);
        }
    }
    else
    {
        std::cout << "\n Error: the variable '" << variable_name << "' to write is not particle data!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
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
void BaseParticles::resizeParticleData<DataType>::
operator()(ParticleData &particle_data, size_t new_size) const
{
    constexpr int type_index = DataTypeIndex<DataType>::value;

    for (size_t i = 0; i != std::get<type_index>(particle_data).size(); ++i)
        std::get<type_index>(particle_data)[i]->resize(new_size, ZeroData<DataType>::value);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addParticleDataWithDefaultValue<DataType>::
operator()(ParticleData &particle_data) const
{
    constexpr int type_index = DataTypeIndex<DataType>::value;

    for (size_t i = 0; i != std::get<type_index>(particle_data).size(); ++i)
        std::get<type_index>(particle_data)[i]->push_back(ZeroData<DataType>::value);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::copyParticleData<DataType>::
operator()(ParticleData &particle_data, size_t index, size_t another_index) const
{
    constexpr int type_index = DataTypeIndex<DataType>::value;

    for (size_t i = 0; i != std::get<type_index>(particle_data).size(); ++i)
        (*std::get<type_index>(particle_data)[i])[index] =
            (*std::get<type_index>(particle_data)[i])[another_index];
}
//=================================================================================================//
template <typename StreamType>
void BaseParticles::writeParticlesToVtk(StreamType &output_stream)
{
    size_t total_real_particles = total_real_particles_;

    // write current/final particle positions first
    output_stream << "   <Points>\n";
    output_stream << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        Vec3d particle_position = upgradeToVec3d(pos_[i]);
        output_stream << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";
    output_stream << "   </Points>\n";

    // write header of particles data
    output_stream << "   <PointData  Vectors=\"vector\">\n";

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

    // compute derived particle variables
    for (auto &derived_variable : derived_variables_)
    {
        derived_variable->exec();
    }

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
void WriteAParticleVariableToXml::
operator()(const std::string &variable_name, StdLargeVec<DataType> &variable) const
{
    SimTK::Xml::element_iterator ele_ite = xml_engine_.root_element_.element_begin();
    for (size_t i = 0; i != total_real_particles_; ++i)
    {
        xml_engine_.setAttributeToElement(ele_ite, variable_name, variable[i]);
        ele_ite++;
    }
}
//=================================================================================================//
template <typename DataType>
void ReadAParticleVariableFromXml::
operator()(const std::string &variable_name, StdLargeVec<DataType> &variable) const
{
    SimTK::Xml::element_iterator ele_ite = xml_engine_.root_element_.element_begin();
    for (size_t i = 0; i != total_real_particles_; ++i)
    {
        xml_engine_.getRequiredAttributeValue(ele_ite, variable_name, variable[i]);
        ele_ite++;
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