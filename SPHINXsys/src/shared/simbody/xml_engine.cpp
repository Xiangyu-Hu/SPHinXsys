/**
 * @file 	xml.cpp
 * @brief 	XML functions are defined here
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "xml_engine.h"

using std::string;

namespace SPH 
{
	//===============================================================//
	XmlEngine::XmlEngine(const std::string &xml_name, const std::string &root_tag)
	{
		xml_name_ = xml_name;
		xmldoc_.setRootTag(root_tag);
  		root_element_ = xmldoc_.getRootElement();
	}
	//===============================================================//
	XmlEngine::~XmlEngine()
	{
		//nothing done here right now
	}
	//===============================================================//
	void XmlEngine::creatXmlElement(const std::string &ele_name)
	{
		element_ = new SimTK::Xml::Element(ele_name);
	}
	//===============================================================//
	void XmlEngine::AddElementToXmlDoc()
	{
		root_element_.insertNodeAfter(root_element_.node_end(), *element_);
	}
	//===============================================================//
	void XmlEngine::AddAttributeToElement(const std::string &attrib_name, const Matd value)
	{
		int num_dim = value.nrow();
		SimTK::Array_<Real> array_(num_dim * num_dim);
		for (int i = 0; i < num_dim; i++)
			for (int j = 0; j < num_dim; j++)
				array_[i * num_dim + j] = value(i, j);
		SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(array_));
		element_->setAttributeValue(attr_.getName(), attr_.getValue());
	}
	//===============================================================//
	void XmlEngine::WriteToXmlFile(const std::string &filefullpath)
	{
		xmldoc_.writeToFile(filefullpath);
	}
	//===============================================================//
	void XmlEngine::LoadXmlFile(const std::string &filefullpath)
	{
		xmldoc_.readFromFile(filefullpath);
		root_element_ = xmldoc_.getRootElement();
	}
	//===============================================================//
	void XmlEngine::ReadFromXmlFile(const std::string &filefullpath)
	{
		//to be written
	}
	//===============================================================//
	std::string XmlEngine::GetRootElementTag()
	{
		return xmldoc_.getRootTag();
	}
	//===============================================================//
	Matd XmlEngine::GetRequiredAttributeMatrixValue(SimTK::Xml::element_iterator &ele_ite_, const std::string &attrib_name)
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
		Matd mat_;
		for (int i = 0; i < num_dim; i++)
			for (int j = 0; j < num_dim; j++)
				mat_(i, j) = array_[i * num_dim + j];
		return mat_;
	}
	//===============================================================//
}