/**
* @file 	case.h
* @brief 	This is the test of using levelset to generate particles relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (3D).
*			Before particle generation, we clean the sharp corner and smooth 0 levelset value, then doing the re-initialization
* @author 	Yongchuan Yu and Xiangyu Hu
*/

#ifndef TEST_3D_PARTICLE_GENERATOR_SHELL_CASE_H
#define TEST_3D_PARTICLE_GENERATOR_SHELL_CASE_H

#include "sphinxsys.h"
#include "particle_generator_lattice_shell.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_airfoil = "./input/plate_50_50_4.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-30, -30, -3);
Vec3d domain_upper_bound(30, 30, 3);
Real dp_0 = 1.5;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

TriangleMeshShape *CreateImportedModelSurface()
{
	Vecd translation(0.0, 0.0, 0.0);
	TriangleMeshShape *geometry_imported_model = new TriangleMeshShape(full_path_to_airfoil, translation, 1.0);

	return geometry_imported_model;
}

class ImportedModel : public SolidBody
{
public:
	ImportedModel(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(), new ShellParticleGeneratorLattice())
	{
		/** Geometry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addTriangleMeshShape(CreateImportedModelSurface(), ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
	}
};

#endif //TEST_3D_PARTICLE_GENERATOR_SHELL_CASE_H