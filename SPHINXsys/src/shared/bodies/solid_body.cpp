/**
 * @file    solid_body.cpp
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "solid_body.h"

#include "base_kernel.h"
#include "sph_system.h"
#include "base_material.h"
#include "solid_particles.h"

namespace SPH {
	//=================================================================================================//
	SolidBody::SolidBody(SPHSystem &system, std::string body_name,
		ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
		: RealBody(system, body_name, particle_adaptation, particle_generator) {}
	//=================================================================================================//
	ThinStructure::ThinStructure(SPHSystem& system, std::string body_name,
		ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
		: SolidBody(system, body_name, particle_adaptation, particle_generator)
	{
		particle_adaptation->getKernel()->reduceOnce();
	}
	//=================================================================================================//
	SolidBodyPartForSimbody
		::SolidBodyPartForSimbody(SPHBody* solid_body, std::string solid_body_part_name)
		: BodyPartByParticle(solid_body, solid_body_part_name)
	{
		solid_particles_ = dynamic_cast<SolidParticles*>(body_->base_particles_);
		Solid* solid = dynamic_cast<Solid*>(body_->base_particles_->base_material_);
		solid_body_density_ = solid->ReferenceDensity();
	}
	//=================================================================================================//
}
