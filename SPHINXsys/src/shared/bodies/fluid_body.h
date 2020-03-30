/**
 * @file    fluid_body.h
 * @brief 	This is the class for bodies used for fluid.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_body.h"
#include "particle_generator_lattice.h"

#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @brief preclaimed class.
	 */
	class SPHSystem;
	/**
	 * @class FluidBody
	 * @brief Declaration of fluid body.
	 */
	class FluidBody : public RealBody
	{
	public:
		explicit FluidBody(SPHSystem &system, string body_name,
			int refinement_level, ParticlesGeneratorOps op);
		virtual ~FluidBody() {};

		/** Build inner configuration. */
		virtual void BuildInnerConfiguration() override;
		/** Build contact configuration. */
		virtual void BuildContactConfiguration() override;
	};
}