/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	base_particles.hpp
 * @brief 	This is the implementation of the template functions in base_particles.h 
 * @author	Xiangyu Hu
 */

#ifndef BASE_PARTICLES_HPP
#define BASE_PARTICLES_HPP

#include "base_particles.h"

namespace SPH
{
    //=================================================================================================//
    template <typename VariableType>
    void BaseParticles::
        registerAVariable(StdLargeVec<VariableType> &variable_addrs,
                          const std::string &variable_name, VariableType initial_value)
    {
        constexpr int type_index = ParticleDataTypeIndex<VariableType>::value;

        if (all_variable_maps_[type_index].find(variable_name) == all_variable_maps_[type_index].end())
        {
            variable_addrs.resize(real_particles_bound_, initial_value);
            std::get<type_index>(all_particle_data_).push_back(&variable_addrs);
            all_variable_maps_[type_index].insert(make_pair(variable_name, std::get<type_index>(all_particle_data_).size() - 1));
        }
        else
        {
            std::cout << "\n Error: the variable '" << variable_name << "' has already been registered!" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }
    //=================================================================================================//
    template <typename VariableType>
    void BaseParticles::
        registerAVariable(StdLargeVec<VariableType> &variable_addrs,
                          const std::string &new_variable_name, const std::string &old_variable_name)
    {
        registerAVariable(variable_addrs, new_variable_name);

        constexpr int type_index = ParticleDataTypeIndex<VariableType>::value;

        if (all_variable_maps_[type_index].find(old_variable_name) != all_variable_maps_[type_index].end())
        {
            StdLargeVec<VariableType> *old_variable =
                std::get<type_index>(all_particle_data_)[all_variable_maps_[type_index][old_variable_name]];
            for (size_t i = 0; i != real_particles_bound_; ++i)
                variable_addrs[i] = (*old_variable)[i];
        }
        else
        {
            std::cout << "\n Error: the old variable '" << old_variable_name << "' is not registered!" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }
    //=================================================================================================//
    template <typename VariableType>
    StdLargeVec<VariableType> *BaseParticles::getVariableByName(std::string variable_name)
    {
        constexpr int type_index = ParticleDataTypeIndex<VariableType>::value;

        if (all_variable_maps_[type_index].find(variable_name) != all_variable_maps_[type_index].end())
            return std::get<type_index>(all_particle_data_)[all_variable_maps_[type_index][variable_name]];

        std::cout << "\n Error: the variable '" << variable_name << "' is not registered!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
        return nullptr;
    }
    //=================================================================================================//
    template <typename VariableType>
    void BaseParticles::
        addAVariableNameToList(ParticleVariableList &variable_name_list, std::string variable_name)
    {
        constexpr int type_index = ParticleDataTypeIndex<VariableType>::value;

        if (all_variable_maps_[type_index].find(variable_name) != all_variable_maps_[type_index].end())
        {
            bool is_to_add = true;
            for (size_t i = 0; i != variable_name_list[type_index].size(); ++i)
            {
                if (variable_name_list[type_index][i].first == variable_name)
                    is_to_add = false;
            }
            if (is_to_add)
            {
                size_t variable_index = all_variable_maps_[type_index][variable_name];
                variable_name_list[type_index].push_back(make_pair(variable_name, variable_index));
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
    template <typename VariableType>
    void BaseParticles::addAVariableToWrite(std::string variable_name)
    {
        addAVariableNameToList<VariableType>(variables_to_write_, variable_name);
    }
    //=================================================================================================//
    template <typename VariableType>
    void BaseParticles::addAVariableToRestart(std::string variable_name)
    {
        addAVariableNameToList<VariableType>(variables_to_restart_, variable_name);
    }
    //=================================================================================================//
    template <typename VariableType>
    void BaseParticles::registerASortableVariable(std::string variable_name)
    {
        constexpr int type_index = ParticleDataTypeIndex<VariableType>::value;

        if (sortable_variable_maps_[type_index].find(variable_name) == sortable_variable_maps_[type_index].end())
        {
            if (all_variable_maps_[type_index].find(variable_name) != all_variable_maps_[type_index].end())
            {
                StdLargeVec<VariableType> *variable =
                    std::get<type_index>(all_particle_data_)[all_variable_maps_[type_index][variable_name]];
                std::get<type_index>(sortable_data_).push_back(variable);
                sortable_variable_maps_[type_index].insert(make_pair(variable_name, std::get<type_index>(sortable_data_).size() - 1));
            }
            else
            {
                std::cout << "\n Error: the variable '" << variable_name << "' for sorting is not registered!" << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cout << "\n Warning: the variable '" << variable_name << "' is already a sortabele variable!" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        }
    }
    //=================================================================================================//
    template <typename VariableType>
    void BaseParticles::addAParticleDataValue<VariableType>::
    operator()(ParticleData &particle_data) const
    {
        constexpr int type_index = ParticleDataTypeIndex<VariableType>::value;

        for (size_t i = 0; i != std::get<type_index>(particle_data).size(); ++i)
            std::get<type_index>(particle_data)[i]->push_back(VariableType(0));
    }
    //=================================================================================================//
    template <typename VariableType>
    void BaseParticles::copyAParticleDataValue<VariableType>::
    operator()(ParticleData &particle_data, size_t this_index, size_t another_index) const
    {
        constexpr int type_index = ParticleDataTypeIndex<VariableType>::value;

        for (size_t i = 0; i != std::get<type_index>(particle_data).size(); ++i)
            (*std::get<type_index>(particle_data)[i])[this_index] =
                (*std::get<type_index>(particle_data)[i])[another_index];
    }
    //=================================================================================================//
    template <typename VariableType>
    void WriteAParticleVariableToXml::
    operator()(std::string &variable_name, StdLargeVec<VariableType> &variable) const
    {
        SimTK::Xml::element_iterator ele_ite = xml_engine_.root_element_.element_begin();
        for (size_t i = 0; i != total_real_particles_; ++i)
        {
            xml_engine_.setAttributeToElement(ele_ite, variable_name, variable[i]);
            ele_ite++;
        }
    }
    //=================================================================================================//
    template <typename VariableType>
    void ReadAParticleVariableFromXml::
    operator()(std::string &variable_name, StdLargeVec<VariableType> &variable) const
    {
        SimTK::Xml::element_iterator ele_ite = xml_engine_.root_element_.element_begin();
        for (size_t i = 0; i != total_real_particles_; ++i)
        {
            xml_engine_.getRequiredAttributeValue(ele_ite, variable_name, variable[i]);
            ele_ite++;
        }
    }
    //=================================================================================================//
}
#endif //BASE_PARTICLES_HPP