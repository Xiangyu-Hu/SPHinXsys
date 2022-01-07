/**
* @file 	3d_self_contact.h
* @brief 	This is the test to check self contact for solid dynamics
* @author 	Xiangyu Hu
*/

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/coil.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real half_width = 55.0;
Real resolution_ref = half_width / 30.0;
Real BW = resolution_ref * 4;
Vec3d domain_lower_bound(-half_width - BW, -half_width - 1.5 * BW, -BW);
Vec3d domain_upper_bound(half_width + BW, half_width + BW, 2.0 * half_width + BW);
// Domain bounds of the system.
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
Real rho0_s = 1.265;
Real poisson = 0.45;
Real Youngs_modulus = 5e4;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	Bodies used in the case.
//----------------------------------------------------------------------
class Coil : public SolidBody
{
public:
	Coil(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd translation(0.0, 0.0, 0.0);
		TriangleMeshShapeSTL triangle_mesh_heart_shape(full_path_to_file, translation, 1.0);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_heart_shape);
	}
};
class StationaryPlate : public SolidBody
{
public:
	StationaryPlate(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize_plate(half_width + BW, 0.5 * BW, half_width + BW);
		Vecd translation_plate(0.0, -half_width - 0.75 * BW, half_width);
		int resolution(20); //SimTK geometric modeling resolution.
		body_shape_.add<TriangleMeshShapeBrick>(halfsize_plate, resolution, translation_plate);
	}
};
