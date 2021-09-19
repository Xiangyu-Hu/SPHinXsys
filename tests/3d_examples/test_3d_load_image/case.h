/**
* @file 	case.h
* @brief 	This is the test of using levelset(distance map) to generate particles relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (3D).
*			Before particle generation, we clean the sharp corner and smooth 0 levelset value, then doing the re-initialization
* @author 	Yijin Mao
*/

#ifndef TEST_3D_LOAD_IMAGE_CASE_H
#define TEST_3D_LOAD_IMAGE_CASE_H

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_image = "./input/sphere.mhd";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-25.0, -25.0, -25.0);
Vec3d domain_upper_bound(25.0, 25.0, 25.0);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 50.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

ImageMeshShape *CreateImportedModelSurface()
{
	Vec3d center(0.0, 0.0, 0.0);
	Vec3d spacings(1.0, 1.0, 1.0);
	ImageMeshShape *geometry_imported_model = \
		new ImageMeshShape(full_path_to_image);

	return geometry_imported_model;
}

class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name,
			new ParticleSpacingByBodyShape(1.15, 1.0, 2),
			new ParticleGeneratorMultiResolution())
	{
		/** Geometry definition. */
		std::unique_ptr<ComplexShapeImageMesh> original_body_shape_mesh(new ComplexShapeImageMesh());
		original_body_shape_mesh->addImageMeshShape(CreateImportedModelSurface(), \
			ShapeBooleanOps::add);
		ComplexShape original_body_shape(original_body_shape_mesh.get());
		body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
	}
};

#endif //TEST_3D_LOAD_IMAGE_CASE_Hs