/**
* @file 	particle_generator_single_resolution.h
* @brief 	This is the test of using levelset to generate particles with single resolution and relax particles.
* @details	We use this case to test the particle generation and relaxation by levelset for a complex geometry (2D).
*			Before particle generation, we clean the sharp corner and smooth 0 levelset value, then doing the reinitialization

* @author 	Yongchuan Yu and Xiangyu Hu
* @version 0.1
*/

#pragma once

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string input_body = "./input/SPHinXsys-2d.dat";

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.3; 				/**< InputBody length right part. */
Real DL1 = 2.3;				/**< InputBody length left part. */
Real DH = 4.5; 				/**< InputBody height. */
Real resolution_ref =(DL + DL1) / 80; 	/**< Reference resolution. */
BoundingBox system_domain_bounds(Vec2d(-DL1, 0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	InputBodyas a solid body
//----------------------------------------------------------------------
class InputBody : public SolidBody
{
public:
	InputBody(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		ComplexShape original_body_shape;
		original_body_shape.addAPolygonFromFile(input_body, ShapeBooleanOps::add);
		body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
	}
};