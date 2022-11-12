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
    /**
     * @class DiscreteVariable
     * @brief template base class for all discrete variables.
     */
    template <typename DataType>
    class DiscreteVariable
    {
        bool isRegistered_;
        size_t index_in_container_;
        const std::sting variable_name_;
        const DataType default_value_;

    public:
        DiscreteVariable(const std::string &variable_name, const DataType &default_value)
            : isRegistered_(false), data_index_in_container_(MaxSize_t),
              variable_name_(variable_name), default_value_(default_value){};
        virtual ~DiscreteVariable(){};

        DataType DefaultValue() { return default_value_; };

        size_t IndexInContainer()
        {
            if (!isRegistered_)
            {
                std::cout << "\n Error: the variable:" << variable_name_ << " is registered!" << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
            return index_in_container_;
        };

        template <class ContainerType, typename InitializationFunction>
        void registerToContainer(ContainerType &container, const InitializationFunction &initialization)
        {
            index_in_container_ = container.registerDiscreteVariable(variable_name_, initialization);
            isRegistered_ = true;
        };
    };
}
#endif // BASE_VARIABLES_H
