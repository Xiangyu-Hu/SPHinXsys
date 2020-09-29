/**
* @file 	sphere.h
* @brief 	This is the case file for the network growing.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
* @version 0.1
*/

#pragma once

#include "sphinxsys.h"
using namespace SPH;

/** Set the file path to the stl file. */
std::string full_path_to_stl = "./input/sphere.stl";
Vec3d domain_lower_bound(-1.0,-1.0, -1.0);
Vec3d domain_upper_bound(1.0, 1.0, 1.0);			
Real dp_0 	= (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
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
	MyPolygonBody(SPHSystem &system, string body_name, int refinement_level, 
		ParticleGenerator* particle_generator)
		: SolidBody(system, body_name, refinement_level, particle_generator)
	{
		ComplexShape original_body_shape;
		original_body_shape.addTriangleMeshShape(CreateCADGeometry(), ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape, 4);
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