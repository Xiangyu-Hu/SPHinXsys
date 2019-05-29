/**
 * @file 	xml_engine.h
 * @brief 	XML class for xml input and output, this is GUI of simbody xml parser.
 * @author	Chi Zhang and Xiangyu Hu.
 * @version	0.1.
 */
#pragma once

#include "base_data_package.h"

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
        SimTK::Xml::Element root_element_;    /**< Root element of docment. */
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
		void AddAttributeToElement(const std::string &attrib_name,const std::string &value);
		/**
  		 * @brief	Adds attribute of type int to an xml element.
  		 * @param[in] 	attrib_name Name of the attribute.
  		 * @param[in] 	value String type value of the attribute.
  		 */
		void AddAttributeToElement(const std::string &attrib_name,const int value);
		/**
  		 * @brief	Adds attribute of type Real to an xml element.
  		 * @param[in] 	attrib_name Name of the attribute.
  		 * @param[in] 	value String type value of the attribute.
  		 */
		void AddAttributeToElement(const std::string &attrib_name,const Real value);
		/**
  		 * @brief	Adds attribute of type vector in 2D to an xml element.
  		 * @param[in] 	attrib_name Name of the attribute.
  		 * @param[in] 	value String type value of the attribute.
  		 */
		void AddAttributeToElement(const std::string &attrib_name,const Vec2d value);
		/**
  		 * @brief	Adds attribute of type vector in 3D to an xml element.
  		 * @param[in] 	attrib_name Name of the attribute.
  		 * @param[in] 	value String type value of the attribute.
  		 */
		void AddAttributeToElement(const std::string &attrib_name,const Vec3d value);
        /**
         * @brief Adds attribute of type matrix to an xml element.
         * @param[in]   attrib_name Name of the attribute.
         * @param[in]   value String type value of the attribute.
        */
        void AddAttributeToElement(const std::string &attrib_name,const Matd value);
		/**
  		 * @brief	Adds attribute of type boolean to an xml element.
  		 * @param[in] 	attrib_name Name of the attribute.
  		 * @param[in] 	value String type value of the attribute.
  		 */
		void AddAttributeToElement(const std::string &attrib_name,const bool value);
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
        /**
         * @brief Get the required int attribute vlaue of an element
         * @param[in] ele input element.
         * @param[in] attrib_name required attribute name.
         * @returns int type value of rquaired attribute.
         */
        int GetRequiredAttributeIntValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name);
        /**
         * @brief Get the required int attribute vlaue of an element
         * @param[in] ele input element.
         * @param[in] attrib_name required attribute name.
         * @returns Real type value of rquaired attribute.
         */
        Real GetRequiredAttributeRealValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name);
        /**
         * @brief Get the required int attribute vlaue of an element
         * @param[in] ele input element.
         * @param[in] attrib_name required attribute name.
         * @returns Vector(in 2D or 3D) type value of rquaired attribute.
         */
        Vecd GetRequiredAttributeVectorValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name);
		/**
         * @brief Get the required int attribute vlaue of an element
         * @param[in] ele input element.
         * @param[in] attrib_name required attribute name.
         * @returns Vector(in 2D or 3D) type value of rquaired attribute.
         */
        Vec3d GetRequiredAttributeVec3dValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name);
		/**
         * @brief Get the required int attribute vlaue of an element
         * @param[in] ele input element.
         * @param[in] attrib_name required attribute name.
         * @returns Vector(in 2D or 3D) type value of rquaired attribute.
         */
        Matd GetRequiredAttributeMatrixValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name);
    };

}