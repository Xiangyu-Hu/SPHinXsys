/*
* @file 	tank_case.h
*/

#ifndef	TANK_CASE_H
#define TANK_CASE_H

#include "sphinxsys.h"
#define PI (3.14159265358979323846)
using namespace SPH;

/* Domain bounds of the system*/
BoundingBox system_domain_bounds(Vec3d(-0.3, -0.3,-0.3), Vec3d(0.3, 0.5,0.3));
Real resolution_ref = 0.01;   /* Initial particle spacing*/

/*
Material properties of the fluid.
*/
Real rho0_f = 1000.0;         /*Fluid density*/
Real rho0_a = 1.0;            /*Air density*/
Real gravity_g = 9.81;        /*Gravity force of fluid*/
Real U_f = 2.0* sqrt(gravity_g * 0.174);	/**< Characteristic velocity. */
Real U_g = 2.0* sqrt(gravity_g * 0.174);  	/**< dispersion velocity in shallow water. */
Real c_f = 10.0 * SMAX(U_g, U_f);	/**< Reference sound speed. */

Real length_scale = 1.0;
Vec3d translation(0, 0.175, 0);
/*
Geometry of the tank, water, air, and sensors.
*/
std::string fuel_tank_outer = "./input/validation_tank_outer_slim.STL";
std::string fuel_tank_inner = "./input/validation_tank_inner.STL";
std::string water_05 = "./input/validation_water.STL";
std::string air_05 = "./input/validation_air.STL";
std::string probe_s1_shape = "./input/ProbeS1.STL";
std::string probe_s2_shape = "./input/ProbeS2.STL";
std::string probe_s3_shape = "./input/ProbeS3.STL";

/*
The Tank.
*/
//class Tank : public ComplexShape
//{
//public:
//	explicit Tank(const std::string &shape_name) :ComplexShape(shape_name)
//	{
//		/** Geometry definition. */
//		
//		add<TriangleMeshShapeSTL>(fuel_tank_outer, translation, length_scale,"OuterWall");
//		subtract<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale,"InnerWall");
//	}
//};
class WallAndStructure : public ComplexShape
{
public:
	explicit WallAndStructure(const std::string& shape_name) :ComplexShape(shape_name)
	{
		/** Geometry definition. */
		add<TriangleMeshShapeSTL>(fuel_tank_inner, translation, length_scale, "InnerWall");
		
	}
};

/*
The Water.
*/
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(water_05, translation, length_scale);
	}
};

/*
The Air.
*/
class AirBlock : public ComplexShape
{
public:
	explicit AirBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(air_05, translation, length_scale);
	}
};

/*
External Excitation
*/
class VariableGravity : public Gravity
{
	Real time_ = 0;
public:
	VariableGravity() : Gravity(Vecd(0.0, -gravity_g, 0.0)) {};
	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		time_= GlobalStaticVariables::physical_time_;
		if (time_ > 0.25)
		{
			global_acceleration_[0] = 4 * PI * PI * 1.63 * 1.63 * 0.0075 * sin(2 * PI * 1.63 * (time_ - 0.25));
			/*global_acceleration_[0] =0;*/
		}
		
		return global_acceleration_;
	}
};

/*
Sensors: S1, S2 and S3;
*/
class ProbeS1 : public ComplexShape
{
public:
	explicit ProbeS1(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s1_shape, translation_probe, length_scale);
	}
};

class ProbeS2 : public ComplexShape
{
public:
	explicit ProbeS2(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_2(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s2_shape, translation_probe_2, length_scale);
	}
};

class ProbeS3 : public ComplexShape
{
public:
	explicit ProbeS3(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vec3d translation_probe_3(0.0,0.0, 0.0);
		add<TriangleMeshShapeSTL>(probe_s3_shape, translation_probe_3, length_scale);
	}
};

class WaterObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit WaterObserverParticleGenerator(SPHBody& sph_body) : ObserverParticleGenerator(sph_body)
	{
		// add observation points
		
	}
};


#endif //TANK_CASE_H