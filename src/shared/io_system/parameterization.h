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
 * @file 	parameterization.h
 * @brief 	This is the base classes for introducing the parameterization
 *			of a class or method.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SPHINXSYS_PARAMETERIZATION_H
#define SPHINXSYS_PARAMETERIZATION_H

#include "base_data_package.h"
#include "xml_parser.h"

namespace SPH
{
/**
 * @class ParameterizationIO
 * @brief Base class for parameter IO.
 */
class ParameterizationIO
{
  public:
    XmlParser xml_parameters_;
    std::string filefullpath_;

    explicit ParameterizationIO(const std::string &input_path);
    ~ParameterizationIO() {};

    void writeProjectParameters();
};
/**
 * @class  BaseParameterization
 * @brief Derived class for parameterization.
 */
template <class BaseClassType>
class BaseParameterization : public BaseClassType
{
  public:
    template <typename... Args>
    explicit BaseParameterization(ParameterizationIO *parameterization_io, Args... args)
        : BaseClassType(std::forward<Args>(args)...),
          xml_parameters_(parameterization_io->xml_parameters_),
          filefullpath_(parameterization_io->filefullpath_){};
    ~BaseParameterization() {};

  protected:
    XmlParser &xml_parameters_;
    std::string filefullpath_;

    template <typename VariableType>
    void getAParameter(const std::string &element_name, const std::string &variable_name, VariableType &variable_addrs)
    {
        auto element = xml_parameters_.findElement(xml_parameters_.first_element_, element_name);
        xml_parameters_.queryAttributeValue(element, variable_name, variable_addrs);
    };

    template <typename VariableType>
    void setAParameter(const std::string &element_name, const std::string &variable_name, VariableType &variable_addrs)
    {
        xml_parameters_.addNewElement(xml_parameters_.first_element_, element_name);
        auto element = xml_parameters_.findElement(xml_parameters_.first_element_, element_name);
        xml_parameters_.setAttributeToElement(element, variable_name, variable_addrs);
    };
};
} // namespace SPH
#endif // SPHINXSYS_PARAMETERIZATION_H
