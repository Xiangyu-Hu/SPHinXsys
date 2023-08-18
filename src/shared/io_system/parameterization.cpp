#include "parameterization.h"

namespace SPH
{
//=============================================================================================//
ParameterizationIO::ParameterizationIO(const std::string &input_path)
    : xml_parameters_("xml_parameters", "parameters"),
      filefullpath_(input_path + "/" + "project_parameters.dat")
{
    xml_parameters_.loadXmlFile(filefullpath_);
}
//=============================================================================================//
void ParameterizationIO::writeProjectParameters()
{
    xml_parameters_.writeToXmlFile(filefullpath_);
}
//=================================================================================================//
} // namespace SPH
