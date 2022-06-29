/**
 * @file 	parameterization.cpp
 * @author	Xiangyu Hu
 */

#include "parameterization.h"

namespace SPH
{
	//=============================================================================================//
	ParameterizationIO::ParameterizationIO(const std::string &input_path)
		: xml_paremeters_("xml_parameters", "parameters"),
		filefullpath_(input_path + "/" + "project_parameters.dat")
	{
		xml_paremeters_.loadXmlFile(filefullpath_);
	}
	//=============================================================================================//
	void ParameterizationIO::writeProjectParameters()
	{
		xml_paremeters_.writeToXmlFile(filefullpath_);
	}
	//=================================================================================================//
}
