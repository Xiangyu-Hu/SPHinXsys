/**
* @file 	airfoil_2d.h
* @brief 	This is the test of using levelset to generate particles in a relax configuration.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (2D).
*			Before particle generation, we clean the sharp corner and smooth zero levelset contour, 
*			which represent the body surface then doing the particle relaxation.
* @author 	Yongchuan Yu and Xiangyu Hu
*/

#ifndef AIRFOIL_2D_H
#define AIRFOIL_2D_H

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string airfoil_flap_front = "./input/airfoil_flap_front.dat";
std::string airfoil_wing = "./input/airfoil_wing.dat";
std::string airfoil_flap_rear = "./input/airfoil_flap_rear.dat";
 //----------------------------------------------------------------------
 //	Basic geometry parameters and numerical setup.
 //----------------------------------------------------------------------
Real DL = 1.25; 				/**< airfoil length rear part. */
Real DL1 = 0.25;				/**< airfoil length front part. */
Real DH = 0.25; 				/**< airfoil height. */
Real resolution_ref = 0.01; 	/**< Reference resolution. */	
BoundingBox system_domain_bounds(Vec2d(-DL1, -DH), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Airfoil	as a solid body
//----------------------------------------------------------------------
class Airfoil : public SolidBody
{
public:
	Airfoil(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name, 
			new ParticleSpacingByBodyShape(1.15, 0, 2), 
			new ParticleGeneratorMultiResolution())
	{
		/** Geometry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addAPolygonFromFile(airfoil_flap_front, ShapeBooleanOps::add);
		original_body_shape.addAPolygonFromFile(airfoil_wing, ShapeBooleanOps::add);
		original_body_shape.addAPolygonFromFile(airfoil_flap_rear, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
	}
};
#endif //AIRFOIL_2D_H