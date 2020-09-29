/**
* @file 	boeing737.cpp
* @brief 	This is the test of using levelset to generate particles relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (3D).
*			Before particle generation, we clean the sharp corner and smooth 0 levelset value, then doing the re-initialization

* @author 	Yongchuan Yu and Xiangyu Hu
* @version 0.1
*/

#pragma once

#include "sphinxsys.h"

using namespace SPH;

/**
 * @brief Basic geometry parameters.
 */
 /** Set the file path to the stl file. */
std::string full_path_to_airfoil = "./input/teapot.stl";

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Vec3d domain_lower_bound(-9.0, -6.0, 0.0);
Vec3d domain_upper_bound(9.0, 6.0, 9.0);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;

TriangleMeshShape *CreateImportedModelSurface()
{
	Vecd translation(0.0, 0.0, 0.0);
	TriangleMeshShape *geometry_imported_model = new TriangleMeshShape(full_path_to_airfoil, translation, 1.0);

	return geometry_imported_model;
}

class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geomerty definition. */
		ComplexShape original_body_shape;
		original_body_shape.addTriangleMeshShape(CreateImportedModelSurface(), ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
	}
};

