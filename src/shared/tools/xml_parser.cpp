#include "xml_parser.h"

namespace SPH
{
//=================================================================================================//
void XmlParser::addNewElement(const std::string &element_name)
{
    xml_doc_->InsertEndChild(xml_doc_->NewElement(element_name.c_str()));
}
//=================================================================================================//
void XmlParser::addNewElement(tinyxml2::XMLElement *element, const std::string &child_name)
{
    element->InsertNewChildElement(child_name.c_str());
}
//=================================================================================================//
size_t XmlParser::Size()
{
    size_t num = 0;
    for (const tinyxml2::XMLElement *child = xml_doc_->FirstChildElement();
         child;
         child = child->NextSiblingElement())
        ++num;

    return num;
}
//=================================================================================================//
size_t XmlParser::Size(tinyxml2::XMLElement *base)
{
    size_t num = 0;
    for (const tinyxml2::XMLElement *child = base->FirstChildElement();
         child;
         child = child->NextSiblingElement())
        ++num;

    return num;
}
//=================================================================================================//
tinyxml2::XMLElement *XmlParser::findElement(const std::string &element_tag)
{
    tinyxml2::XMLElement *child_element = xml_doc_->FirstChildElement(element_tag.c_str());
    return child_element;
}
//=================================================================================================//
tinyxml2::XMLElement *XmlParser::findElement(tinyxml2::XMLElement *base, const std::string &element_tag)
{
    tinyxml2::XMLElement *child_element = base->FirstChildElement(element_tag.c_str());
    return child_element;
}
//=================================================================================================//
void XmlParser::resize(const size_t input_size, const std::string name)
{
    size_t total_elements = XmlParser::Size();

    if (total_elements <= input_size)
    {
        for (size_t i = total_elements; i != input_size; ++i)
            XmlParser::addNewElement(name);
    }
    else
    {
        auto child = xml_doc_->FirstChildElement();
        for (size_t i = input_size; i != total_elements; ++i)
        {
            auto delete_child = child->NextSiblingElement();
            xml_doc_->DeleteChild(delete_child);
        }
    }
}
//=================================================================================================//
void XmlParser::resize(tinyxml2::XMLElement *element, const size_t input_size, const std::string name)
{
    size_t total_elements = XmlParser::Size(element);

    if (total_elements <= input_size)
    {
        for (size_t i = total_elements; i != input_size; ++i)
            XmlParser::addNewElement(element, name);
    }
    else
    {
        auto child = element->FirstChildElement();
        for (size_t i = input_size; i != total_elements; ++i)
        {
            auto delete_child = child->NextSiblingElement();
            element->DeleteChild(delete_child);
        }
    }
}
//=================================================================================================//
} // namespace SPH
