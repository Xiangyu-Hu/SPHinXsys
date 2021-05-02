/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file    solid_body.h
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
  */

#pragma once

#include "base_body.h"

namespace SPH {
	/**
	 * @brief Preclaimed class.
	 */
	class SPHSystem;
	class SolidParticles;
	/**
	 * @class SolidBody
	 * @brief Declaration of solidbody which is used for Solid BCs and derived from RealBody.
	 */
	class SolidBody : public RealBody
	{
	public:
		SolidBody(SPHSystem &system, std::string body_name,
			ParticleAdaptation* particle_adaptation = new ParticleAdaptation(1.15),
			ParticleGenerator* particle_generator = new ParticleGeneratorLattice());
		virtual ~SolidBody() {};
		virtual SolidBody* ThisObjectPtr()  override {return this;};
	};

	/**
	 * @class ThinStructure
	 * @brief Declaration of thin structure solidbody.
	 */
	class ThinStructure : public SolidBody
	{
	public:
		ThinStructure(SPHSystem& system, std::string body_name,
			ParticleAdaptation* particle_adaptation = new ParticleAdaptation(1.15),
			ParticleGenerator* particle_generator = new ParticleGeneratorLattice());
		virtual ~ThinStructure() {};
		virtual ThinStructure* ThisObjectPtr() override {return this;};
	};

	/**
	 * @class SolidBodyPartForSimbody
	 * @brief A SolidBodyPart for coupling with Simbody.
	 * The mass, origin, and unit inertial matrix are computed.
	 * Note: In Simbody, all spatial vectors are three dimensional.
	 */
	class SolidBodyPartForSimbody : public BodyPartByParticle
	{
	public:
		Vec3d initial_mass_center_;
		SimTK::MassProperties* body_part_mass_properties_;
		
		SolidBodyPartForSimbody(SPHBody* body, std::string solid_body_part_name);
		virtual~SolidBodyPartForSimbody() {};
	protected:
		Real solid_body_density_;
		SolidParticles* solid_particles_;

		virtual void tagBodyPart() override;
	};	
}
