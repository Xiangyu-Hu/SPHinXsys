/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	base_variables.h
 * @brief 	Here gives classes for the base variables used in simulation.
 * @details These variables are those discretized in spaces and time.
 * @author	Xiangyu Hu
 */

#ifndef BASE_VARIABLES_H
#define BASE_VARIABLES_H

#include "base_data_package.h"

namespace SPH
{
    template <typename DataType>
    class DiscreteVariable;
    const bool sharedVariable = true;
    typedef DataContainerAddressAssemble<DiscreteVariable> DiscreteVariableAssemble;

    /**
     * @class DiscreteVariable
     * @brief template base class for all discrete variables.
     */
    template <typename DataType>
    class DiscreteVariable
    {
    public:
        DiscreteVariable(DiscreteVariableAssemble &variable_assemble,
                         const std::string &name, bool is_shared = !sharedVariable)
            : name_(name), index_in_container_(initializeIndex(variable_assemble, is_shared)){};
        virtual ~DiscreteVariable(){};
        size_t IndexInContainer() const { return index_in_container_; };
        std::string VariableName() const { return name_; };

    private:
        const std::string name_;
        size_t index_in_container_;

        size_t initializeIndex(DiscreteVariableAssemble &variable_assemble, bool is_shared)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            auto &variable_container = std::get<type_index>(variable_assemble);
            size_t determined_index = determineIndex(variable_container);

            if (determined_index == variable_container.size()) // determined a new index
            {
                variable_container.push_back(this);
            }
            else if (!is_shared)
            {
                std::cout << "\n Error: the variable: " << name_ << " is already used!" << std::endl;
                std::cout << "\n Please check if " << name_ << " is a sharable variable." << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }

            return determined_index;
        };

        template <typename VariableContainer>
        size_t determineIndex(const VariableContainer &variable_container)
        {
            size_t i = 0;
            while (i != variable_container.size())
            {
                if (variable_container[i]->VariableName() == name_)
                {
                    return i;
                }
                ++i;
            }
            return variable_container.size();
        }
    };
}
#endif // BASE_VARIABLES_H
