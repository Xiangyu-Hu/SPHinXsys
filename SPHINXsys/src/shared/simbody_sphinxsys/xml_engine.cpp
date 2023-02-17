#include "xml_engine.h"

namespace SPH
{
	//=================================================================================================//
	XmlEngine::XmlEngine(const std::string &xml_name, const std::string &root_tag) : xml_name_(xml_name)
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
		const std::string& attrib_name, const Real& value)
	{
		SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(value));
		ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
	}
	//=================================================================================================//
	void XmlEngine::setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite, 
		const std::string& attrib_name, const int& value)
	{
		SimTK::Xml::Attribute attr_(attrib_name, SimTK::String(value));
		ele_ite->setAttributeValue(attr_.getName(), attr_.getValue());
	}
	//=================================================================================================//
	void XmlEngine::setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite, 
		const std::string& attrib_name, const Vec2d& value)
	{
		setAttributeToVectorElement<2>(ele_ite, attrib_name, value);
	}
	//=================================================================================================//
	void XmlEngine::setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite,
		const std::string& attrib_name, const Vec3d& value)
	{
		setAttributeToVectorElement<3>(ele_ite, attrib_name, value);
	}
	void XmlEngine::setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite, 
		const std::string& attrib_name, const Mat2d& value)
	{
		setAttributeToMatrixElement<2>(ele_ite, attrib_name, value);
	}
	//=================================================================================================//
	void XmlEngine::setAttributeToElement(const SimTK::Xml::element_iterator& ele_ite, 
		const std::string& attrib_name, const Mat3d& value)
	{
		setAttributeToMatrixElement<3>(ele_ite, attrib_name, value);
	}
	//=================================================================================================//
	void XmlEngine::getRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, 
		const std::string& attrib_name, Real& value)
	{
		std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
		value = SimTK::convertStringTo<Real>(value_in_string);
	}
	//=================================================================================================//
	void XmlEngine::getRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, 
		const std::string& attrib_name, int& value)
	{
		std::string value_in_string = ele_ite_->getRequiredAttributeValue(attrib_name);
		value = SimTK::convertStringTo<int>(value_in_string);
	}
	//=================================================================================================//
	void XmlEngine::getRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, 
		const std::string& attrib_name, Vec2d& value)
	{
		getVectorAttributeValue<2>(ele_ite_, attrib_name, value);
	}
	//=================================================================================================//
	void XmlEngine::getRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, 
		const std::string& attrib_name, Vec3d& value)
	{
		getVectorAttributeValue<3>(ele_ite_, attrib_name, value);
	}
	//=================================================================================================//
	void XmlEngine::getRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, 
		const std::string& attrib_name, Mat2d& value)
	{
		getMatrixAttributeValue<2>(ele_ite_, attrib_name, value);
	}
	//=================================================================================================//
	void XmlEngine::getRequiredAttributeValue(SimTK::Xml::element_iterator& ele_ite_, 
		const std::string& attrib_name, Mat3d& value)
	{
		getMatrixAttributeValue<3>(ele_ite_, attrib_name, value);
	}
	//=================================================================================================//
	void XmlEngine::writeToXmlFile(const std::string &filefullpath)
	{
		xmldoc_.writeToFile(filefullpath);
	}
	//=================================================================================================//
	void XmlEngine::loadXmlFile(const std::string &filefullpath)
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
	std::string XmlEngine::getElementTag(SimTK::Xml::Element &element)
	{
		return element.getElementTag();
	}
	//=================================================================================================//
	void XmlEngine::resizeXmlDocForParticles(size_t input_size)
	{
		size_t total_elements = std::distance(root_element_.element_begin(),
											  root_element_.element_end());

		if (total_elements <= input_size)
		{
			for (size_t i = total_elements; i != input_size; ++i)
				addElementToXmlDoc("particle");
		}
		else
		{
			std::cout << "\n Error: XML Engine allows increase date size only!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
	size_t XmlEngine::SizeOfXmlDoc()
	{
		return std::distance(root_element_.element_begin(), root_element_.element_end());
	}
	//=================================================================================================//
	SimTK::Xml::Element XmlEngine::getChildElement(const std::string &tag)
	{
		return root_element_.getOptionalElement(tag);
	}
	//=================================================================================================//
}
