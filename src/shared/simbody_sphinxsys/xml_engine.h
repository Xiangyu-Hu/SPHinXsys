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
 * @file 	xml_engine.h
 * @brief 	XML class for xml input and output, this is GUI of simbody xml parser.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef XML_ENGINE_SIMBODY_H
#define XML_ENGINE_SIMBODY_H

#include "SimTKcommon/internal/Xml.h"
#include "base_data_package.h"
#include "simbody_middle.h"
#include "sph_data_containers.h"
#include "type_wrapper.h"

#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

namespace SPH
{
class XmlEngine
{
  protected:
    std::string xml_name_;        /**< xml name. */
    SimTK::Xml::Document xmldoc_; /**< the xml document. */

  public:
    /** Constructor for XML output.  */
    XmlEngine(const std::string &xml_name, const std::string &root_tag);
    /** Default destructor. */
    virtual ~XmlEngine(){};

    SimTK::Xml::Element root_element_; /**< Root element of document. */

    /**Add existing element to root_element of Xml Doc. */
    void addElementToXmlDoc(const std::string &element_name);

    /**Add child element to a given element. */
    void addChildToElement(SimTK::Xml::Element &father_element, const std::string &child_name);

    //----------------------------------------------------------------------
    //	Add an attribute of type string to an xml element.
    //----------------------------------------------------------------------
    template <typename T>
    void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name, const T &value)
    {
        SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(value));
        ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
    };

    template <int DIMENSION, auto... Rest>
    void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name,
                               const Eigen::Matrix<Real, DIMENSION, 1, Rest...> &value)
    {
        SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(EigenToSimTK(value)));
        ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
    };
    template <int DIMENSION, auto... Rest>
    void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name,
                               const Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
    {
        SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(EigenToSimTK(value)));
        ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
    };

    //----------------------------------------------------------------------
    //	Get the required attribute value of an element.
    //----------------------------------------------------------------------
    template <typename T>
    void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name, T &value)
    {
        std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
        value = SimTK::convertStringTo<T>(value_in_string);
    };

    template <int DIMENSION, auto... Rest>
    void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name,
                                   Eigen::Matrix<Real, DIMENSION, 1, Rest...> &value)
    {
        std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
        value = SimTKToEigen(SimTK::convertStringTo<SimTK::Vec<DIMENSION>>(value_in_string));
    };

    template <int DIMENSION, auto... Rest>
    void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name,
                                   Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
    {
        std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
        value = SimTKToEigen(SimTK::convertStringTo<SimTK::Mat<DIMENSION, DIMENSION>>(value_in_string));
    };

    /** Write to XML file */
    void writeToXmlFile(const std::string &filefullpath);
    /**  Load XML file using XML parser. */
    void loadXmlFile(const std::string &filefullpath);
    /** Get the Tag of root element as a string */
    std::string getRootElementTag();
    /** Get the Tag of a element as a string */
    std::string getElementTag(SimTK::Xml::Element &element);
    /** resize of Xml doc */
    void resizeXmlDocForParticles(size_t input_size);
    /** Get the size of Xml doc */
    size_t SizeOfXmlDoc();
    /** Get a reference to a child element */
    SimTK::Xml::Element getChildElement(const std::string &tag);
};

/**
 * @class XMLMemoryIO
 * @brief The base class defines xml memory data operation.
 */
class XmlMemoryIO
{
  public:
    XmlMemoryIO(){};
    virtual ~XmlMemoryIO(){};

    template <typename T>
    void writeDataToXmlMemory(XmlEngine &xml_engine, SimTK::Xml::Element &element, const BiVector<T> &quantity,
                              int snapshot_, int observation_, const std::string &quantity_name, StdVec<std::string> &element_tag)
    {
        for (int snapshot_index = 0; snapshot_index != snapshot_; ++snapshot_index)
        {
            std::string element_name = element_tag[snapshot_index];
            xml_engine.addChildToElement(element, element_name);
            for (int observation_index = 0; observation_index != observation_; ++observation_index)
            {
                SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
                std::string attribute_name_ = quantity_name + "_" + std::to_string(observation_index);
                xml_engine.setAttributeToElement(ele_ite, attribute_name_, quantity[snapshot_index][observation_index]);
            }
        }
    };

    template <typename T>
    void writeDataToXmlMemory(XmlEngine &xml_engine, SimTK::Xml::Element &element,
                              std::string element_name, int observation_index, const T &quantity, const std::string &quantity_name)
    {
        SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
        std::string attribute_name_ = quantity_name + "_" + std::to_string(observation_index);
        xml_engine.setAttributeToElement(ele_ite, attribute_name_, quantity);
    };

    template <typename T>
    void readDataFromXmlMemory(XmlEngine &xml_engine, SimTK::Xml::Element &element,
                               int observation_index, BiVector<T> &result_container, const std::string &quantity_name)
    {
        int snapshot_index = 0;
        SimTK::Xml::element_iterator ele_ite = element.element_begin();
        for (; ele_ite != element.element_end(); ++ele_ite)
        {
            std::string attribute_name_ = quantity_name + "_" + std::to_string(observation_index);
            xml_engine.getRequiredAttributeValue(ele_ite, attribute_name_, result_container[snapshot_index][observation_index]);
            snapshot_index++;
        }
    };

    void readTagFromXmlMemory(SimTK::Xml::Element &element, StdVec<std::string> &element_tag)
    {
        size_t snapshot_index = 0;
        SimTK::Xml::element_iterator ele_ite = element.element_begin();
        for (; ele_ite != element.element_end(); ++ele_ite)
        {
            element_tag[snapshot_index] = ele_ite->getElementTag();
            snapshot_index++;
        }
    };
};
} // namespace SPH

#endif // XML_ENGINE_SIMBODY_H