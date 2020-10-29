/**
 * @file    solid_body.cpp
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "solid_body.h"
#include "sph_system.h"
#include "base_material.h"
#include "solid_particles.h"

namespace SPH {
	//=================================================================================================//
	SolidBody::SolidBody(SPHSystem &system, string body_name,
		int refinement_level, ParticleGenerator* particle_generator)
		: RealBody(system, body_name, refinement_level, 1.05, particle_generator)
	{
	
	}
	//=================================================================================================//
	SolidBodyPartForSimbody
		::SolidBodyPartForSimbody(SPHBody* solid_body, string solid_body_part_name)
		: BodyPartByParticle(solid_body, solid_body_part_name)
	{
		solid_particles_ = dynamic_cast<SolidParticles*>(body_->base_particles_);
		Solid* solid = dynamic_cast<Solid*>(body_->base_particles_->base_material_);
		solid_body_density_ = solid->ReferenceDensity();
	}
	//=================================================================================================//
}
