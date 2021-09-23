/**
* @file 	plastic_case.h
* @brief 	This is the case setup for plastic taylor bar.cpp.
* @author 	xiaojing tang Chi Zhang and Xiangyu Hu
*/
#pragma once

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Global geometry parameters.
//----------------------------------------------------------------------
Real PL = 0.0032; 		/**< X-direction domain. */
Real PW = 0.0324;		/**< Z-direction domain. */
Real particle_spacing_ref = PL / 10.0;
/** YOU can try PW = 0.2 and particle_spacing_ref = PH / 10.0 to see an interesting test. */
Real SL = particle_spacing_ref * 4.0; 	/**< Length of the holder is one layer particle. */
Real inner_circle_radius = PL;

Vec3d domain_lower_bound(-4.0 * PL, -4.0 * PL, -SL);
Vec3d domain_upper_bound(4.0 * PL, 4.0 * PL, 2.0 * PW);
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);
//----------------------------------------------------------------------
//	Material properties and global parameters
//----------------------------------------------------------------------
Real rho0_s = 8930.0; 			/**< Reference density. */
Real poisson = 0.35; 			/**< Poisson ratio. */
Real Youngs_modulus = 1.17e11;
Real yield_stress = 0.4e9;
Real hardening_modulus = 0.1e9;

TriangleMeshShape* CreateColumn()
{
	Vecd translation_column(0.0, 0.0, 0.6 * PW);
	TriangleMeshShape* geometry = new TriangleMeshShape(SimTK::UnitVec3(0, 0, 1.0), inner_circle_radius,
		0.5 * PW, resolution, translation_column);
	return geometry;
}

/** Define the body constrain geometry. */
TriangleMeshShape* CreateHolder()
{
	Vecd halfsize_holder(3.0 * PL, 3.0 * PL, 0.5 * SL);
	Vecd translation_holder(0.0, 0.0, -0.5 * SL);
	TriangleMeshShape* geometry_holder = new TriangleMeshShape(halfsize_holder,
		resolution, translation_holder);
	return geometry_holder;
}

class Wall : public SolidBody
{
public:
	Wall(SPHSystem& system, string body_name) :
		SolidBody(system, body_name, new ParticleAdaptation(1.15, 1.0))
	{
		ComplexShapeTriangleMesh *mesh = new ComplexShapeTriangleMesh();
		body_shape_ = new ComplexShape(mesh);
		mesh->addTriangleMeshShape(CreateHolder(), ShapeBooleanOps::add);
	}
};


/** Define the body. */
class Column : public SolidBody
{
public:
	Column(SPHSystem& system, std::string body_name) : 
		SolidBody(system, body_name, new ParticleAdaptation(1.3, 1.0))
	{
		ComplexShapeTriangleMesh *mesh = new ComplexShapeTriangleMesh();
		ComplexShape original_body_shape(mesh);
		mesh->addTriangleMeshShape(CreateColumn(), ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape);
	}
};

class WallMaterial : public LinearElasticSolid
{
public:
	WallMaterial() :LinearElasticSolid()
	{
		rho0_ = rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};

/**
 * Setup body material property
 */
class PlasticColumnMaterial : public HardeningPlasticSolid
{
public:
	PlasticColumnMaterial() : HardeningPlasticSolid()
	{
		rho0_ = rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		yield_stress_ = yield_stress;
		hardening_modulus_ = hardening_modulus;
		assignDerivedMaterialParameters();
	}
};
/**
 * application dependent initial condition
 */

class InitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	InitialCondition(SolidBody* body)
		: solid_dynamics::ElasticDynamicsInitialCondition(body) {};
protected:
	void Update(size_t index_i, Real dt) override
	{
			vel_n_[index_i][2] = -227.0;
	}
};

//define an observer body
class MyObserver : public FictitiousBody
{
public:
	MyObserver(SPHSystem& system, string body_name)
		: FictitiousBody(system, body_name)
	{
		body_input_points_volumes_.push_back(make_pair(Vecd(0.0, 0.0, PW), 0.0));
		body_input_points_volumes_.push_back(make_pair(Vecd(PL, 0.0, 0.0), 0.0));
	}
};
