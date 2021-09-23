/**
* @file 	sphere.h
* @brief 	This is the case file for the network growing.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Chi Zhang and Xiangyu Hu
*/

#ifndef TEST_3D_PARTICLE_GENERATION_SPHERE_H
#define TEST_3D_PARTICLE_GENERATION_SPHERE_H

#include "sphinxsys.h"
using namespace SPH;

/** Set the file path to the stl file. */
std::string full_path_to_stl = "./input/sphere.stl";
Vec3d domain_lower_bound(-1.0,-1.0, -1.0);
Vec3d domain_upper_bound(1.0, 1.0, 1.0);			
Real dp_0 	= (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

/** Define the geometry. */
TriangleMeshShape* CreateCADGeometry()
{
	Vecd translation(-1.0, -1.0, -1.0);
	TriangleMeshShape *geometry = new TriangleMeshShape(full_path_to_stl, translation, 0.025);

	return geometry;
}
/** Define the my heart body. */
class MyPolygonBody : public SolidBody
{
public:
	MyPolygonBody(SPHSystem &system, string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(1.15, 1.0),
			new ParticleGeneratorNetwork(Vecd(-1.0, 0.0, 0.0), Vecd(-0.964, 0.0, 0.266), 15, 5.0))
	{
		ComplexShapeTriangleMesh *mesh = new ComplexShapeTriangleMesh();
		ComplexShape original_body_shape(mesh);
		mesh->addTriangleMeshShape(CreateCADGeometry(), ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
/** the material for the newwork. */
class BodyMaterial : public LinearElasticSolid
{
public:
	BodyMaterial() : LinearElasticSolid()
	{
		rho0_ = 1.0;
		E0_ = 1.0;
		nu_ = 0.0;

		assignDerivedMaterialParameters();
	}
};
#endif //TEST_3D_PARTICLE_GENERATION_SPHERE_H