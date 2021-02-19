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
 * @file 	xml_engine.h
 * @brief 	XML class for xml input and output, this is GUI of simbody xml parser.
 * @author	Chi Zhang and Xiangyu Hu.
 * @version	0.1.
 */
#pragma once
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_package.h"

#include "SimTKcommon/internal/Xml.h"
#include "SimTKcommon/internal/String.h"

#include <iostream>
#include <string>
#include <cstdio>

#include <fstream>
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

namespace SPH {
	/**
	 * @class 	XmlEngine
	 * @brief 	XmlEngine class, in which SIMBody XML parse is used.
	 */

	class XmlEngine
	{
		std::string xml_name_;				  /**< xml name. */
	public:
        /**
         * @brief   Defaut constructor.
         */
        XmlEngine(){};
		/**
	 	 * @brief 	Constructor for XML output.
	 	 * @param[in]	xml_name 	The name the Xml object.
	 	 * @param[in]	root_tag 	The name of the root tag.
	 	 */
  		XmlEngine(const std::string &xml_name, const std::string &root_tag);
		/**
	 	 * @brief 	Defaut distructor.
	 	 */
  		~XmlEngine();
        SimTK::Xml::Document xmldoc_;         /**< the xml document. */
        SimTK::Xml::Element root_element_;    /**< Root element of document. */
        SimTK::Xml::Element* element_;        /**< pointer to element. */
        //SimTK::Xml::element_iterator ele_ite_; /**< element iterator. */
  		/**
	 	 * @brief 	Creat an Xml Element.
	 	 * @param[in]	ele_name Name of the element.
	 	 */
  		void CreatXmlElement(const std::string &ele_name);
  		/**
  		 * @brief	Add existing element to root_element of Xml Doc.
  		 */
  		void AddElementToXmlDoc();
		/**
		 * @brief	Adds attribute of type string to an xml element.
		 * @param[in] 	attrib_name  Name of the attribute.
		 * @param[in] 	value_ 	String type value of the attribute.
		 */
		template<class T>
		void AddAttributeToElement(const std::string& attrib_name, const T& value) {
			SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(value));
			element_->setAttributeValue(attr_.getName(), attr_.getValue());
		};
		/**
		* @brief Adds attribute of type matrix to an xml element.
		* @param[in]   attrib_name Name of the attribute.
		* @param[in]   value String type value of the attribute.
		*/
		void AddAttributeToElement(const std::string &attrib_name, const Matd value);
		/**
		  * @brief Get the required int attribute vlaue of an element
		  * @param[in] ele input element.
		  * @param[in] attrib_name required attribute name.
		  * @returns Vector(in 2D or 3D) type value of rquaired attribute.
		  */
		template<class T>
		T GetRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, const std::string& attrib_name) {
			std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
			return SimTK::convertStringTo<T>(value_in_string);
		};
		/**
		* @brief Get the required int attribute valaue of an element
		* @param[in] ele input element.
		* @param[in] attrib_name required attribute name.
		* @returns Vector(in 2D or 3D) type value of required attribute.
		*/
		Matd GetRequiredAttributeMatrixValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name);
		/**
  		 * @brief	Write to XML file
  		 * @param[in] 	filefullpath  Full path to writting file.
  		 */
  		void WriteToXmlFile(const std::string &filefullpath);
         /**
         * @brief   Load XML file using XML parser.
         * @param[in]   filefullpath  Full path to writting file.
         */
        void LoadXmlFile(const std::string &filefullpath);
         /**
         * @brief   Read from XML file
         * @param[in]   filefullpath  Full path to writting file.
         */
        void ReadFromXmlFile(const std::string &filefullpath);
        /**
         * @ Get the Tag of root element as a string
         * @returns Tag of root element.
         */
        std::string GetRootElementTag();
    };

}