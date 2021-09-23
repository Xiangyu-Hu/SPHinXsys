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
 */

#ifndef XML_ENGINE_SIMBODY_H
#define XML_ENGINE_SIMBODY_H


#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_package.h"
#include "sph_data_containers.h"

#include "SimTKcommon.h"
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

namespace SPH 
{
	class XmlEngine
	{
	protected:
		std::string xml_name_;			/**< xml name. */
		SimTK::Xml::Document xmldoc_;	/**< the xml document. */

	public:
		/** Constructor for XML output.  */
		XmlEngine(const std::string& xml_name, const std::string& root_tag);
		/** Defaut distructor. */
		virtual ~XmlEngine() {};

		SimTK::Xml::Element root_element_;	/**< Root element of document. */

		/**Add existing element to root_element of Xml Doc. */
		void addElementToXmlDoc(const std::string& element_name);

		/**Add child element to a given element. */
		void addChildToElement(SimTK::Xml::Element& father_element, const std::string& child_name);

		/** Add an attribute of type string to an xml element.  */
		template<class T>
		void setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite, const std::string& attrib_name, const T& value)
		{
			SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(value));
			ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
		};
		/** Adds attribute of type matrix to an xml element. */
		void setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite, const std::string& attrib_name, const Matd& value);

		/** Get the required attribute value of an element */
		template<class T>
		void getRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, const std::string& attrib_name, T& value)
		{
			std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
			value = SimTK::convertStringTo<T>(value_in_string);
		};
		/** Get the required int attribute valaue of an element */
		void getRequiredAttributeMatrixValue(SimTK::Xml::element_iterator& ele_ite_, const std::string& attrib_name, Matd& value);

		/** Write to XML file */
		void writeToXmlFile(const std::string& filefullpath);
		/**  Load XML file using XML parser. */
		void loadXmlFile(const std::string& filefullpath);
		/** Get the Tag of root element as a string */
		std::string getRootElementTag();
		/** Get the Tag of a element as a string */
		std::string getElementTag(SimTK::Xml::Element& element);
		/** resize of Xml doc */
		void resizeXmlDocForParticles(size_t input_size);
		/** Get the size of Xml doc */
		size_t SizeOfXmlDoc();
		/** Get a reference to a child element */
		SimTK::Xml::Element getChildElement(const std::string& tag);
	};
}

#endif //XML_ENGINE_SIMBODY_H