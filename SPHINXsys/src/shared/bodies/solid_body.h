/**
 * @file    solid_body.h
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#pragma once

#include "base_body.h"

using namespace std;

namespace SPH {
	/**
	 * @brief Preclaimed class.
	 */
	class SPHSystem;
	/**
	 * @class SolidBody
	 * @brief Declaration of solidbody which is used for Solid BCs and derived from RealBody.
	 */
	class SolidBody : public RealBody
	{
	public:
		SolidBody(SPHSystem &system, string body_name, 
			int refinement_level, ParticlesGeneratorOps op);
		virtual ~SolidBody() {};

		/** Set up the reaction model, if essential */
		/** The pointer to derived class object. */
		virtual SolidBody* PointToThisObject() override { return this; };
	};
}