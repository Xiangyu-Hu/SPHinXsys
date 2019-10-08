/**
 * @file    fluid_body.h
 * @brief 	This is the class for bodies used for fluid.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_body.h"
#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @brief preclaimed class.
	 */
	class SPHSystem;
	class Fluid;
	class FluidParticles;
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

	/**
	  * @class FluidBodyPart
	  * @brief An auxillariy class for fluidBody to
	  * indicate a part of the body fixed with location.
	  */
	class FluidBodyPart : public EulerianBodyPart
	{
	protected:
		FluidBody *fluid_body_;
		/** The geometric region indicating the body part. */
		Region fluid_body_part_region_;

		virtual void TagBodyPartCells() override;
	public:
		FluidBodyPart(FluidBody *fluid_body, string fluid_body_part_name);
		virtual~FluidBodyPart() {};
	};
}