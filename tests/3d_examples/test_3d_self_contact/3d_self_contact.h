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
std::string full_path_to_airfoil = "./input/coil.stl";
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
Real rho_0 = 1.265;
Real poisson = 0.45;
Real Youngs_modulus = 5e4;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	Geometries used in the case.
//----------------------------------------------------------------------
TriangleMeshShape *createImportedModelSurface()
{
	Vecd translation(0.0, 0.0, 0.0);
	TriangleMeshShape *geometry_imported_model = new TriangleMeshShape(full_path_to_airfoil, translation, 1.0);

	return geometry_imported_model;
}
int resolution(20); //SimTK geometric modeling resolution.
TriangleMeshShape* createStationaryPlate()
{
	Vecd halfsize_plate(half_width + BW, 0.5 * BW, half_width + BW);
	Vecd translation_plate(0.0, -half_width - 0.75 * BW, half_width);
	TriangleMeshShape* geometry_plate = new TriangleMeshShape(halfsize_plate,
		resolution, translation_plate);

	return geometry_plate;
}
//----------------------------------------------------------------------
//	Bodies used in the case.
//----------------------------------------------------------------------
class Coil : public SolidBody
{
public:
	Coil(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		ComplexShapeTriangleMesh *mesh = new ComplexShapeTriangleMesh();
		ComplexShape original_body_shape(mesh);
		mesh->addTriangleMeshShape(createImportedModelSurface(), ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);

	}
};
class StationaryPlate : public SolidBody
{
public:
	StationaryPlate(SPHSystem& system, std::string body_name)
		: SolidBody(system, body_name)
	{
		mesh_.reset(new ComplexShapeTriangleMesh());
		body_shape_ = new ComplexShape(mesh_.get());
		mesh_->addTriangleMeshShape(createStationaryPlate(), ShapeBooleanOps::add);
	}
private:
	std::unique_ptr<ComplexShapeTriangleMesh> mesh_;
};
//----------------------------------------------------------------------
//	Materials used in the case.
//----------------------------------------------------------------------
class CoilMaterial : public NeoHookeanSolid
{
public:
	CoilMaterial() : NeoHookeanSolid()
	{
		rho0_ = rho_0;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};
class StationaryPlateMaterial : public LinearElasticSolid
{
public:
	StationaryPlateMaterial() : LinearElasticSolid()
	{
		rho0_ = rho_0;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};
