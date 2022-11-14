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
 * @brief 	This is the base class variables used in simulation.
 * @details These variables are those discretized in spaces and time.
 * @author	Xiangyu Hu and Chi Zhang
 */

#ifndef BASE_VARIABLES_H
#define BASE_VARIABLES_H

#include "base_data_package.h"

namespace SPH
{
    template <typename DataType>
    class DiscreteVariable;

    typedef GeneralDataAssemble<DiscreteVariable> DiscreteVariableAssemble;

    /**
     * @class DiscreteVariable
     * @brief template base class for all discrete variables.
     */
    template <typename DataType>
    class DiscreteVariable
    {
        const std::string name_;
        const DataType initialize_value_;
        size_t index_in_container_;

    public:
        DiscreteVariable(DiscreteVariableAssemble &variable_assemble,
                         const std::string &name, const DataType &value)
            : name_(name), initialize_value_(value),
              index_in_container_(MaxSize_t)
        {
            addTo(variable_assemble);
        };
        virtual ~DiscreteVariable(){};
        DataType InitializeValue() const { return initialize_value_; };
        size_t IndexInContainer() const { return index_in_container_; };

    protected:
        void addTo(DiscreteVariableAssemble &variable_assemble)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            auto &variable_container = std::get<type_index>(variable_assemble);
            index_in_container_ = variable_container.size();
            variable_container.push_back(this);
        };
    };
}
#endif // BASE_VARIABLES_H
