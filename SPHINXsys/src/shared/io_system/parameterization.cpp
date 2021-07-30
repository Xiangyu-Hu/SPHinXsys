/**
 * @file 	parameterization.cpp
 * @author	Xiangyu Hu
 */

#include "parameterization.h"

namespace SPH
{
	//=============================================================================================//
	ParameterizationIO::ParameterizationIO(In_Output& in_output) :
		xml_paremeters_("xml_parameters", "parameters")
	{
		filefullpath_ = in_output.input_folder_ + "/" + "project_parameters.dat";
		xml_paremeters_.loadXmlFile(filefullpath_);
	}
	//=============================================================================================//
	void ParameterizationIO::writeProjectParameters()
	{
		xml_paremeters_.writeToXmlFile(filefullpath_);
	}
	//=================================================================================================//
}
