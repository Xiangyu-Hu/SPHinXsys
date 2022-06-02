/**
 * @file 	xml.hpp
 * @brief 	XML functions are defined here
 * @author	Bo Zhang and Xiangyu Hu
 */

#ifndef XML_ENGINE_SIMBODY_HPP
#define XML_ENGINE_SIMBODY_HPP

#include "xml_engine.h"

namespace SPH
{
	//=================================================================================================//
	template <typename T>
	void XmlMemoryIO::writeDataToXmlMemory(XmlEngine &xmlengine, SimTK::Xml::Element &element, const DoubleVec<T> &quantity,
		int snapshot_, int observation_, const std::string &quantity_name, StdVec<std::string> &element_tag)
	{
		for (int snapshot_index = 0; snapshot_index != snapshot_; ++snapshot_index)
		{
			std::string element_name = element_tag[snapshot_index];
			xmlengine.addChildToElement(element, element_name);
			for (int observation_index = 0; observation_index != observation_; ++observation_index)
			{
				SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
				std::string attribute_name_ = quantity_name + "_" + std::to_string(observation_index);
				xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity[snapshot_index][observation_index]);
			}
		}
	};
	//=================================================================================================//
	template <typename T>
	void XmlMemoryIO::writeDataToXmlMemory(XmlEngine &xmlengine, SimTK::Xml::Element &element,
		std::string element_name, int observation_index, const T &quantity, const std::string &quantity_name)
	{
		SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
		std::string attribute_name_ = quantity_name + "_" + std::to_string(observation_index);
		xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity);
	};
	//=================================================================================================//
	template <typename T>
	void XmlMemoryIO::readDataFromXmlMemory(XmlEngine &xmlengine, SimTK::Xml::Element &element,
		int observation_index, DoubleVec<T> &result_container, const std::string &quantity_name)
	{
		int snapshot_index = 0;
		SimTK::Xml::element_iterator ele_ite = element.element_begin();
		for (; ele_ite != element.element_end(); ++ele_ite)
		{
			std::string attribute_name_ = quantity_name + "_" + std::to_string(observation_index);
			xmlengine.getRequiredAttributeValue<T>(ele_ite, attribute_name_, result_container[snapshot_index][observation_index]);
			snapshot_index++;
		}
	};
	//=================================================================================================// 
	void XmlMemoryIO::readDataFromXmlMemory(XmlEngine &xmlengine, SimTK::Xml::Element &element,
		size_t observation_index, DoubleVec<Matd> &result_container, const std::string &quantity_name)
	{
		size_t snapshot_index = 0;
		SimTK::Xml::element_iterator ele_ite = element.element_begin();
		for (; ele_ite != element.element_end(); ++ele_ite)
		{
			std::string attribute_name_ = quantity_name + "_" + std::to_string(observation_index);
			xmlengine.getRequiredAttributeMatrixValue(ele_ite, attribute_name_, result_container[snapshot_index][observation_index]);
			snapshot_index++;
		}
	};
	//=================================================================================================//
	void XmlMemoryIO::readTagFromXmlMemory(SimTK::Xml::Element &element, StdVec<std::string> &element_tag)
	{
		size_t snapshot_index = 0;
		SimTK::Xml::element_iterator ele_ite = element.element_begin();
		for (; ele_ite != element.element_end(); ++ele_ite)
		{
			element_tag[snapshot_index] = ele_ite->getElementTag();
			snapshot_index++;
		}
	};
	//=================================================================================================//
	
	//=================================================================================================//
}
#endif //XML_ENGINE_SIMBODY_HPP
