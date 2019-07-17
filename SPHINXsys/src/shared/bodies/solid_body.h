/**
 * @file    solid_body.h
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#pragma once

#include "base_body.h"
#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @brief Preclaimed class.
	 */
	class SPHSystem;
	class Solid;
	class SolidParticles;

	/**
	 * @class SolidBody
	 * @brief Declaration of solidbody which is used for Solid BCs and derived from RealBody.
	 */
	class SolidBody : public RealBody
	{
	public:
		SolidBody(SPHSystem &system, string body_name, Solid &solid_material,
			SolidParticles &solid_particles, int refinement_level, ParticlesGeneratorOps op);
		virtual ~SolidBody() {};

		/** material of this body. */
		Solid &solid_material_;
		/** Particles in this body. */
		SolidParticles &solid_particles_;	

		/** Build inner configuration. */
		virtual void BuildInnerConfiguration() override;
		/** Build contact configuration. */
		virtual void BuildContactConfiguration() override;
		/** The pointer to derived class object. */
		virtual SolidBody* PointToThisObject() override { return this; };
	};

	/**
	 * @class SolidBodyPart
	 * @brief A auxillariy class for SoildBody to
	 * indicate a part of the body moving together with particles.
	 */
	class SolidBodyPart : public LagrangianBodyPart
	{
	protected:
		SolidBody *solid_body_;
		Region soild_body_part_region_;

		virtual void TagBodyPartParticles() override;
	public:
		SolidBodyPart(SolidBody *solid_body, string soild_body_part_name);
		virtual~SolidBodyPart() {};
	};

	/**
	 * @class SolidBodyPartForSimbody
	 * @brief A SolidBodyPart for coupling with Simbody.
	 * The mass, origin, and unit inertial matrix are computed
	 * Note: Simbody spatial vectors are three dimensional
	 */
	class SolidBodyPartForSimbody : public SolidBodyPart
	{
	protected:
		Real solid_body_density_;


		virtual void TagBodyPartParticles() override;
	public:
		SolidBodyPartForSimbody(SolidBody *solid_body,
			string soild_body_part_name, Real solid_body_density);
		virtual~SolidBodyPartForSimbody() {};

		SimTK::MassProperties *body_part_mass_properties_;
		Vec3 initial_mass_center_;
	};
}