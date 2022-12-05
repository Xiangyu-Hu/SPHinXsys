/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	xml_engine.h
 * @brief 	XML class for xml input and output, this is GUI of simbody xml parser.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef XML_ENGINE_SIMBODY_H
#define XML_ENGINE_SIMBODY_H

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "array.h"

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
#include <filesystem>
namespace fs = std::filesystem;
#endif

namespace SPH
{
	class XmlEngine
	{
	protected:
		std::string xml_name_;		  /**< xml name. */
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

		/** Add an attribute of type string to an xml element.  */
		template <class T>
		void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name, const T &value);

		void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name, const int &value)
		{
			SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(value));
			ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
		};
		void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name, const Real &value)
		{
			SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(value));
			ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
		};
		void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name, const Veci &value)
		{
			SimTK::Array_<int, int> array_(Dimensions);

			for (int i = 0; i < Dimensions; ++i)
				array_[i] = value(i);

			SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(array_));
			ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
		};
		void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name, const Vecd &value)
		{
			SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(EigenToSimTK(value)));
			ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
		};
		void setAttributeToElement(const SimTK::Xml::element_iterator &ele_ite, const std::string &attrib_name, const Matd &value)
		{
			SimTK::Array_<Real, int> array_(Dimensions * Dimensions);

			for (int i = 0; i < Dimensions; i++)
				for (int j = 0; j < Dimensions; j++)
					array_[i * Dimensions + j] = value(i, j);

			SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(array_));
			ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
		};
		/** Get the required attribute value of an element */
		template <class T>
		void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name, T &value);
		/** Get the required int attribute value of an element */
		void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name, int &value)
		{
			std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
			value = SimTK::convertStringTo<int>(value_in_string);
		};
		/** Get the required real attribute value of an element */
		void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name, Real &value)
		{
			std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
			value = SimTK::convertStringTo<Real>(value_in_string);
		};
		/** Get the required int vector attribute value of an element */
		void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name, Veci &value)
		{
			std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
			SimTK::Array_<int, int> array_;
			array_ = SimTK::convertStringTo<SimTK::Array_<int, int>>(value_in_string);

			if (array_.size() != Dimensions)
			{
				std::cout << "\n Error: the dimension of data in XML is not valid" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (int i = 0; i < Dimensions; i++)
				value[i] = array_[i];
		};
		/** Get the required real vector attribute value of an element */
		void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name, Vecd &value)
		{
			std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
			value = SimTKToEigen(SimTK::convertStringTo<SimTKVecd>(value_in_string));
		};
		/** Get the required matrix attribute value of an element */
		void getRequiredAttributeValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name, Matd &value)
		{
			std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
			SimTK::Array_<Real, int> array_;
			array_ = SimTK::convertStringTo<SimTK::Array_<Real, int>>(value_in_string);

			if (array_.size() != Dimensions * Dimensions)
			{
				std::cout << "\n Error: the dimension of data in XML is not valid" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			for (int i = 0; i < Dimensions; i++)
				for (int j = 0; j < Dimensions; j++)
					value(i, j) = array_[i * Dimensions + j];
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
		void writeDataToXmlMemory(XmlEngine &xml_engine, SimTK::Xml::Element &element, const DoubleVec<T> &quantity,
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
								   int observation_index, DoubleVec<T> &result_container, const std::string &quantity_name)
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
}

#endif // XML_ENGINE_SIMBODY_H