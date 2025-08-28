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
 * @file    xml_parser.h
 * @brief   XML parser based on tinyxml2 (https://github.com/leethomason/tinyxml2)
 * @author	Chi Zhang
 */
#pragma once

#include "tinyxml2.h"

#include "base_data_package.h"
#include "sphinxsys_containers.h"

#include <cassert>
#include <charconv>
#include <cstdio>
#include <iostream>
#include <string>

#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
#define assertm(exp, msg) assert(((void)msg, exp))

namespace SPH
{
// Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]");
// Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
// Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]");
// Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
// Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

template <typename DataType>
inline std::string DataToString(const DataType &value)
{
    std::ostringstream out;
    out.precision(15);
    out << std::fixed << value;
    return std::move(out).str();
    // return std::to_string(value);
}

template <int DIMENSION, auto... Rest>
inline std::string DataToString(const Eigen::Matrix<Real, DIMENSION, Rest...> &value)
{
    std::stringstream ss;
    ss << value.format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", ""));
    return ss.str();
}

template <int DIMENSION, auto... Rest>
inline std::string DataToString(const Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
{
    std::stringstream ss;
    ss << value.format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", ""));
    return ss.str();
}

template <typename DataType>
inline void StringToData(std::string &value_str, DataType &value)
{
    std::istringstream(value_str) >> value;
}

template <int DIMENSION, auto... Rest>
inline void StringToData(std::string &value_str, Eigen::Matrix<Real, DIMENSION, 1, Rest...> &value)
{
    std::vector<Real> temp;
    temp.resize(DIMENSION);
    std::istringstream value_stream(value_str);

    int i = 0;
    for (std::string line; std::getline(value_stream, line, ',');)
    {
        temp[i] = std::stod(line);
        i++;
    }
    assert(i == DIMENSION);

    for (int j = 0; j != DIMENSION; j++)
        value[j] = temp[j];
}

template <int DIMENSION, auto... Rest>
inline void StringToData(std::string &value_str, Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
{
    std::vector<Real> temp;
    temp.resize(DIMENSION * DIMENSION);
    std::istringstream value_stream(value_str);

    int i = 0;
    for (std::string line; std::getline(value_stream, line, ',');)
    {
        temp[i] = std::stod(line);
        i++;
    }

    assert(i == DIMENSION * DIMENSION);

    for (int j = 0; j != DIMENSION; j++)
        for (int k = 0; k != DIMENSION; k++)
            value(j, k) = temp[j * DIMENSION + k];
}

/**
 * @class XmlParser
 * @note The XmlParser represents a wrapper of tinyxml2.
 */

class XmlParser
{
  protected:
    std::string xml_name_;           /**< xml name. */
    tinyxml2::XMLDocument *xml_doc_; /**< the xml document. */

  public:
    /** Constructor and destructor.  */
    XmlParser(const std::string &xml_name)
    {
        xml_name_ = xml_name;
        xml_doc_ = new tinyxml2::XMLDocument();
    };
    XmlParser(const std::string &xml_name, const std::string &first_element_name)
        : XmlParser(xml_name)
    {
        addNewElement(first_element_name);
        first_element_ = xml_doc_->FirstChildElement();
    };
    virtual ~XmlParser() { delete xml_doc_; };

    /**
     * Functions & Parameters
     */
    /** First element of the xml doc, also noted as root elemeng in tinyxml-1. */
    tinyxml2::XMLElement *first_element_;

    /** Write to XML file */
    void writeToXmlFile(const std::string &filefullpath)
    {
        xml_doc_->SaveFile(filefullpath.c_str());
        if (xml_doc_->Error())
        {
            std::cout << "Error in save file: " << xml_doc_->ErrorStr() << std::endl;
            exit(EXIT_FAILURE);
        }
    };

    /**  Load XML file using XML parser. */
    void loadXmlFile(const std::string &filefullpath)
    {
        xml_doc_->LoadFile(filefullpath.c_str());
        if (xml_doc_->Error())
        {
            std::cout << "Error in load file: " << xml_doc_->ErrorStr() << std::endl;
            exit(EXIT_FAILURE);
        }
        first_element_ = xml_doc_->FirstChildElement();
    };

    /**  Parse a xml string. */
    void parseXmlStream(const std::string &xml_stream)
    {
        xml_doc_->Parse(xml_stream.c_str());
        first_element_ = xml_doc_->FirstChildElement();
    };

    /**  Parse a xml string. */
    std::string getXmlErrorName(const tinyxml2::XMLError &xml_error)
    {
        return std::string(xml_doc_->ErrorIDToName(xml_error));
    };

    /** Get the Tag of root element as a string */
    std::string getFirstElementTag()
    {
        const char *first_element_name = first_element_->Name();
        return std::string(first_element_name);
    };

    /** Get the Tag of a element as a string */
    std::string getElementTag(tinyxml2::XMLElement *element)
    {
        const char *element_name = element->Name();
        return std::string(element_name);
    };

    /** Get a reference to a child element */
    tinyxml2::XMLElement *getChildElement(const std::string &tag)
    {
        return findElement(tag);
    };

    /**Add new element to Xml Doc. */
    inline void addNewElement(const std::string &element_name);

    /**Add child element to a given element. */
    inline void addNewElement(tinyxml2::XMLElement *father_element, const std::string &child_name);

    /** Get the size of Xml doc */
    inline size_t Size();

    /** Get the size of Xml Element */
    inline size_t Size(tinyxml2::XMLElement *base);

    /** Find optional element in root element */
    inline tinyxml2::XMLElement *findElement(const std::string &element_tag);

    /** Find optional element in optional element */
    inline tinyxml2::XMLElement *findElement(tinyxml2::XMLElement *base, const std::string &element_tag);

    /** resize of Xml doc */
    inline void resize(const size_t input_size, const std::string name);

    /** resize of an element */
    inline void resize(tinyxml2::XMLElement *element, const size_t input_size, const std::string name);

    //----------------------------------------------------------------------
    //	Add an attribute of type string to an xml element.
    //----------------------------------------------------------------------
    template <typename T>
    void setAttributeToElement(tinyxml2::XMLElement *base_ele, const std::string &attrib_name, const T &value)
    {
        base_ele->SetAttribute(attrib_name.c_str(), DataToString(value).c_str());
    };

    template <int DIMENSION, auto... Rest>
    void setAttributeToElement(tinyxml2::XMLElement *base_ele, const std::string &attrib_name,
                               const Eigen::Matrix<Real, DIMENSION, 1, Rest...> &value)
    {
        base_ele->SetAttribute(attrib_name.c_str(), DataToString(value).c_str());
    };

    template <int DIMENSION, auto... Rest>
    void setAttributeToElement(tinyxml2::XMLElement *base_ele, const std::string &attrib_name,
                               const Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
    {
        base_ele->SetAttribute(attrib_name.c_str(), DataToString(value).c_str());
    };

    //----------------------------------------------------------------------
    //	Get the required attribute value of an element.
    //----------------------------------------------------------------------
    template <typename T>
    void queryAttributeValue(tinyxml2::XMLElement *base_ele, const std::string &attrib_name, T &value)
    {
        const char *value_char = 0;
        base_ele->QueryAttribute(attrib_name.c_str(), &value_char);
        std::string value_str(value_char);

        StringToData(value_str, value);
    };

    template <int DIMENSION, auto... Rest>
    void queryAttributeValue(tinyxml2::XMLElement *base_ele, const std::string &attrib_name,
                             Eigen::Matrix<Real, DIMENSION, 1, Rest...> &value)
    {
        const char *value_char = 0;
        base_ele->QueryAttribute(attrib_name.c_str(), &value_char);
        std::string value_str(value_char);

        StringToData(value_str, value);
    };

    template <int DIMENSION, auto... Rest>
    void queryAttributeValue(tinyxml2::XMLElement *base_ele, const std::string &attrib_name,
                             Eigen::Matrix<Real, DIMENSION, DIMENSION, Rest...> &value)
    {
        const char *value_char = 0;
        base_ele->QueryAttribute(attrib_name.c_str(), &value_char);
        std::string value_str(value_char);

        StringToData(value_str, value);
    };
};

/**Add new element to Xml Doc. */
inline void XmlParser::addNewElement(const std::string &element_name)
{
    xml_doc_->InsertEndChild(xml_doc_->NewElement(element_name.c_str()));
}

/**Add child element to a given element. */
inline void XmlParser::addNewElement(tinyxml2::XMLElement *element, const std::string &child_name)
{
    element->InsertNewChildElement(child_name.c_str());
}

/** Get the size of Xml doc */
inline size_t XmlParser::Size()
{
    size_t num = 0;
    for (const tinyxml2::XMLElement *child = xml_doc_->FirstChildElement();
         child;
         child = child->NextSiblingElement())
        ++num;

    return num;
}

/** Get the size of element child */
inline size_t XmlParser::Size(tinyxml2::XMLElement *base)
{
    size_t num = 0;
    for (const tinyxml2::XMLElement *child = base->FirstChildElement();
         child;
         child = child->NextSiblingElement())
        ++num;

    return num;
}

/** Find optional element in root element */
inline tinyxml2::XMLElement *XmlParser::findElement(const std::string &element_tag)
{
    tinyxml2::XMLElement *child_element = xml_doc_->FirstChildElement(element_tag.c_str());
    return child_element;
}

/** Find optional element in optional element */
inline tinyxml2::XMLElement *XmlParser::findElement(tinyxml2::XMLElement *base, const std::string &element_tag)
{
    tinyxml2::XMLElement *child_element = base->FirstChildElement(element_tag.c_str());
    return child_element;
}

/** resize of Xml doc */
inline void XmlParser::resize(const size_t input_size, const std::string name)
{
    size_t total_elements = XmlParser::Size();

    if (total_elements <= input_size)
    {
        for (size_t i = total_elements; i != input_size; ++i)
            XmlParser::addNewElement(name);
    }
    else
    {
        std::cout << "\n Error: XML Parser allows increase date size only!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}

/** resize of an element */
inline void XmlParser::resize(tinyxml2::XMLElement *element, const size_t input_size, const std::string name)
{
    size_t total_elements = XmlParser::Size(element);

    if (total_elements <= input_size)
    {
        for (size_t i = total_elements; i != input_size; ++i)
            XmlParser::addNewElement(element, name);
    }
    else
    {
        std::cout << "\n Error: XML Parser allows increase date size only!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
} // namespace SPH
