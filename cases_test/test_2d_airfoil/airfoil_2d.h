/**
* @file 	airfoil_2d.cpp
* @brief 	This is the test of using levelset to generate particles relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (2D).
*			Before particle generation, we clean the sharp corner and smooth 0 levelset value, then doing the reinitialization

* @author 	Yongchuan Yu and Xiangyu Hu
* @version 0.1
*/

#pragma once

#include "sphinxsys.h"

using namespace SPH;

/** Set the file path to the dat file. */
std::string airfoil_flap_front = "./input/airfoil_flap_front.dat";
std::string airfoil_wing = "./input/airfoil_wing.dat";
std::string airfoil_flap_rear = "./input/airfoil_flap_rear.dat";

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 5.0; 					/**< Channel length rear part. */
Real DL1 = 3.0;                /**< Channel length front part. */
Real DH = 3.0; 						/**< Channel height. */
Real particle_spacing_ref = 0.005; 	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4.0; 		/**< Boundary width, determined by specific layer of boundary particles. */

/** Airfoil	as a solid body */
class Airfoil : public SolidBody
{
public:
	Airfoil(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geometry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addAPolygonFromFile(airfoil_flap_front, ShapeBooleanOps::add);
		original_body_shape.addAPolygonFromFile(airfoil_wing, ShapeBooleanOps::add);
		original_body_shape.addAPolygonFromFile(airfoil_flap_rear, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
	}
};
