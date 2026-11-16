#include "xml_parser.h"

namespace SPH
{
//=================================================================================================//
void XmlParser::addNewElement(const std::string &element_name)
{
    xml_doc_->InsertEndChild(xml_doc_->NewElement(element_name.c_str()));
}
//=================================================================================================//
void XmlParser::addNewElement(TinyXMLElement *element, const std::string &child_name)
{
    element->InsertNewChildElement(child_name.c_str());
}
//=================================================================================================//
size_t XmlParser::Size()
{
    size_t num = 0;
    for (const TinyXMLElement *child = xml_doc_->FirstChildElement();
         child;
         child = child->NextSiblingElement())
        ++num;

    return num;
}
//=================================================================================================//
size_t XmlParser::Size(TinyXMLElement *base)
{
    size_t num = 0;
    for (const TinyXMLElement *child = base->FirstChildElement();
         child;
         child = child->NextSiblingElement())
        ++num;

    return num;
}
//=================================================================================================//
TinyXMLElement *XmlParser::findElement(const std::string &element_tag)
{
    TinyXMLElement *child_element = xml_doc_->FirstChildElement(element_tag.c_str());
    return child_element;
}
//=================================================================================================//
TinyXMLElement *XmlParser::findElement(TinyXMLElement *base, const std::string &element_tag)
{
    TinyXMLElement *child_element = base->FirstChildElement(element_tag.c_str());
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
void XmlParser::resize(TinyXMLElement *element, const size_t input_size, const std::string name)
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
