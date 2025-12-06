#include <gtest/gtest.h>
#include "sphinxsys.h"
#include "xml_parser.h"

using namespace SPH;


TEST(test_xml, test_xmler)
{
	XmlParser xml_parser( "xml_parser", "element1");
	xml_parser.addNewElement( "element2" );
	xml_parser.addNewElement( xml_parser.first_element_, "element11" );
	xml_parser.addNewElement( xml_parser.findElement( "element2" ), "element22" );

	std::string root_name = xml_parser.getFirstElementTag();

	xml_parser.writeToXmlFile( "output/test_xml_parser.xml" );

	EXPECT_EQ(root_name, "element1");
}

TEST(test_xml, test_writer)
{
	XmlParser xml_parser( "xml_parser");
	xml_parser.loadXmlFile( "input/test.xml" );
	
	std::string root_name = xml_parser.getFirstElementTag();

	xml_parser.writeToXmlFile( "output/simple.xml" );

	EXPECT_EQ(root_name, "document");
}

TEST(test_xml, test_passer)
{
	static const std::string xml =
		"<information>"
		"</information>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	
	std::string root_name = xml_parser.getFirstElementTag();

	EXPECT_EQ(root_name, "information");
}

TEST(test_xml, test_geter)
{
	XmlParser xml_parser( "xml_parser");
	xml_parser.loadXmlFile( "input/test.xml" );
	/** find child element in the first element. */
	std::string ele_name = xml_parser.getElementTag( xml_parser.findElement( xml_parser.first_element_, "English" ) );

	EXPECT_EQ(ele_name, "English");
}

TEST(test_xml, test_doc_size)
{
	static const std::string xml =
		"<information>"
		"</information>"
		"<name>"
		"</name>"
		"<address>"
		"</address>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	
	size_t size_of_doc = xml_parser.Size();

	EXPECT_EQ(size_of_doc, 3);
}

TEST(test_xml, test_element_size)
{
	static const std::string xml =
		"<information>"
		"</information>"
		"<name>"
		"</name>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );

	xml_parser.addNewElement( xml_parser.findElement( "name" ), "first_name" );
	xml_parser.addNewElement( xml_parser.findElement( "name" ), "middle_name" );
	xml_parser.addNewElement( xml_parser.findElement( "name" ), "last_name" );

	size_t size_of_element = xml_parser.Size( xml_parser.findElement( "name" ) );

	EXPECT_EQ(size_of_element, 3);
}

TEST(test_xml, test_doc_resizer)
{
	static const std::string xml =
		"<information>"
		"</information>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	xml_parser.resize(10, "information");

	xml_parser.writeToXmlFile( "output/particle.xml" );
	
	size_t size_of_doc = xml_parser.Size();
	EXPECT_EQ(size_of_doc, 10);
}

TEST(test_xml, test_element_resizer)
{
	static const std::string xml =
		"<information>"
		"</information>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	xml_parser.resize(xml_parser.first_element_, 100, "name");

	xml_parser.writeToXmlFile( "output/information.xml" );
	
	size_t size_of_doc = xml_parser.Size(xml_parser.first_element_);
	EXPECT_EQ(size_of_doc, 100);
}

TEST(test_xml, test_scalar_attributer)
{
	static const std::string xml =
		"<root>"
		"</root>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	xml_parser.resize(xml_parser.first_element_, 100, "particle");

	size_t id = 0;
	for( auto child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		xml_parser.setAttributeToElement( child, "Index", id );
		++id;
	}

	xml_parser.writeToXmlFile( "output/particle_id.xml" );
	
	size_t size_of_doc = xml_parser.Size(xml_parser.first_element_);
	EXPECT_EQ(size_of_doc, 100);
}

TEST(test_xml, test_vector_attributer)
{
	XmlParser xml_parser( "xml_parser", "particles");
	xml_parser.resize(xml_parser.first_element_, 100, "particle");

	for( tinyxml2::XMLElement *child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		Vec2d postion(-3.8, 1.6);
		xml_parser.setAttributeToElement( child, "Position", postion );
	}

	xml_parser.writeToXmlFile( "output/particle_pos.xml" );
	
	size_t size_of_doc = xml_parser.Size(xml_parser.first_element_);
	EXPECT_EQ(size_of_doc, 100);
}

TEST(test_xml, test_matrix_attributer)
{
	static const std::string xml =
		"<root>"
		"</root>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	xml_parser.resize(xml_parser.first_element_, 100, "particle");

	for( tinyxml2::XMLElement *child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		Matd stress = Matd::Random();
		xml_parser.setAttributeToElement( child, "Sigma", stress );
	}

	xml_parser.writeToXmlFile( "output/particle_stress.xml" );
	
	size_t size_of_doc = xml_parser.Size(xml_parser.first_element_);
	EXPECT_EQ(size_of_doc, 100);
}


TEST(test_xml, test_scalar_query)
{
	static const std::string xml =
		"<root>"
		"</root>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	xml_parser.resize(xml_parser.first_element_, 20, "particle");

	for( tinyxml2::XMLElement *child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		xml_parser.setAttributeToElement( child, "Density", 5.0 );
	}

	for( tinyxml2::XMLElement *child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		Real density = 0;
		xml_parser.queryAttributeValue( child, "Density", density );
		EXPECT_EQ( density, 5.0 );
	}
}

TEST(test_xml, test_vector_query)
{
	XmlParser xml_parser( "xml_parser");
	xml_parser.loadXmlFile( "output/particle_pos.xml" );
	
	std::string root_name = xml_parser.getFirstElementTag();


	Vec2d postion(-3.8, 1.6);

	for( tinyxml2::XMLElement *child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		Vec2d pos_query = Vec2d::Zero();
		xml_parser.queryAttributeValue( child, "Position", pos_query );
		EXPECT_EQ( postion, pos_query );
	}
}

TEST(test_xml, test_matrix_query)
{
	static const std::string xml =
		"<root>"
		"</root>";

	XmlParser xml_parser("xml_parser");
	xml_parser.parseXmlStream( xml );
	xml_parser.resize(xml_parser.first_element_, 1, "particle");

	Matd stress = Matd::Identity();
	stress(1,1) = 9.0;


	for( tinyxml2::XMLElement *child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		xml_parser.setAttributeToElement( child, "Position", stress );
	}

	for( tinyxml2::XMLElement *child = xml_parser.first_element_->FirstChildElement(); 
		 child; 
		 child = child->NextSiblingElement())
	{
		Matd str_query = Matd::Zero();
		xml_parser.queryAttributeValue( child, "Position", str_query );
		EXPECT_EQ( stress, str_query );
	}
}


int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
