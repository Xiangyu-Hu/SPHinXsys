/**
 * @file 	xml.cpp
 * @brief 	XML functions are defined here
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "xml_engine.h"

namespace SPH
{
	//=================================================================================================//
	XmlEngine::XmlEngine(const std::string& xml_name, const std::string& root_tag) :
		xml_name_(xml_name)
	{
		xmldoc_.setRootTag(root_tag);
		root_element_ = xmldoc_.getRootElement();
	}
	//=================================================================================================//
	void XmlEngine::addElementToXmlDoc(const std::string &element_name)
	{
		root_element_.insertNodeAfter(root_element_.node_end(), SimTK::Xml::Element(element_name));
	}
	//=================================================================================================//
	void XmlEngine::addChildToElement(SimTK::Xml::Element &father_element,
									  const std::string &child_name)
	{
		father_element.insertNodeAfter(father_element.node_end(), SimTK::Xml::Element(child_name));
	}
	//=================================================================================================//
	void XmlEngine::setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite,
		const std::string& attrib_name, const Matd& value)
	{
		int num_dim = value.nrow();
		SimTK::Array_<Real> array_(num_dim * num_dim);
		for (int i = 0; i < num_dim; i++)
			for (int j = 0; j < num_dim; j++)
				array_[i * num_dim + j] = value(i, j);
		SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(array_));
		ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
	}
	//=================================================================================================//
	void XmlEngine::getRequiredAttributeMatrixValue(SimTK::Xml::element_iterator& ele_ite_, const std::string& attrib_name, Matd& value)
	{
		std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
		SimTK::Array_<Real> array_;
		array_ = SimTK::convertStringTo<SimTK::Array_<float>>(value_in_string);
		int num_dim_2 = array_.size();
		int num_dim;
		if (num_dim_2 == 4) {
			num_dim = 2;
		}
		else if (num_dim_2 == 9) {
			num_dim = 3;
		}
		else
		{
			std::cout << "\n Error: the input dimension of deformation tensor:" << num_dim_2 << " is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		for (int i = 0; i < num_dim; i++)
			for (int j = 0; j < num_dim; j++)
				value(i, j) = array_[i * num_dim + j];
	}
	//=================================================================================================//
	void XmlEngine::writeToXmlFile(const std::string& filefullpath)
	{
		xmldoc_.writeToFile(filefullpath);
	}
	//=================================================================================================//
	void XmlEngine::loadXmlFile(const std::string& filefullpath)
	{
		xmldoc_.readFromFile(filefullpath);
		root_element_ = xmldoc_.getRootElement();
	}
	//=================================================================================================//
	std::string XmlEngine::getRootElementTag()
	{
		return xmldoc_.getRootTag();
	}
	//=================================================================================================//
	std::string XmlEngine::getElementTag(SimTK::Xml::Element& element)
	{
		return element.getElementTag();
	}
	//=================================================================================================//
	void  XmlEngine::resizeXmlDocForParticles(size_t input_size)
	{
		size_t total_elements =  std::distance(root_element_.element_begin(),
			root_element_.element_end());

		if (total_elements <= input_size) 
		{
			for (size_t i = total_elements; i != input_size; ++i) addElementToXmlDoc("particle");
		}
		else
		{
			std::cout << "\n Error: XML Engine allows increase date size only!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	};
	//=================================================================================================//
	size_t XmlEngine::SizeOfXmlDoc()
	{
		return std::distance(root_element_.element_begin(), root_element_.element_end());
	}
	//=================================================================================================//
	SimTK::Xml::Element XmlEngine::getChildElement(const std::string& tag)
	{
		return root_element_.getOptionalElement(tag);
	};
	//=================================================================================================//
}
