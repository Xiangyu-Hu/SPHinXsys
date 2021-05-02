/**
* @file 	sphere.h
* @brief 	This is the case file for the network growing.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
*/

#pragma once

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
	MyPolygonBody(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name, new ParticleAdaptation(1.15, 0),
			new ParticleGeneratorNetwork(Vecd(-1.0, 0.0, 0.0), Vecd(-0.964, 0.0, 0.266)))
	{
		ComplexShape original_body_shape;
		original_body_shape.addTriangleMeshShape(CreateCADGeometry(), ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};
/** the material for the newwork. */
class BodyMaterial : public LinearElasticSolid
{
public:
	BodyMaterial() : LinearElasticSolid()
	{
		rho_0_ = 1.0;
		E_0_ = 1.0;
		nu_ = 0.0;

		assignDerivedMaterialParameters();
	}
};