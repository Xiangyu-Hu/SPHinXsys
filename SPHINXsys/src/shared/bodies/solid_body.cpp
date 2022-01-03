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

namespace SPH
{
	//=================================================================================================//
	SolidBody::SolidBody(SPHSystem &system, const std::string &body_name,
						 SharedPtr<SPHAdaptation> sph_adaptation_ptr)
		: RealBody(system, body_name, sph_adaptation_ptr)
	{
		sph_system_.addASolidBody(this);
	}
	//=================================================================================================//
	ThinStructure::ThinStructure(SPHSystem &system, const std::string &body_name,
								 SharedPtr<SPHAdaptation> sph_adaptation_ptr)
		: SolidBody(system, body_name, sph_adaptation_ptr)
	{
		sph_adaptation_ptr->getKernel()->reduceOnce();
	}
	//=================================================================================================//
	SolidBodyPartForSimbody::
		SolidBodyPartForSimbody(SPHBody &body, const std::string &body_part_name, Shape &shape)
		: BodyRegionByParticle(body, body_part_name, shape),
		  solid_particles_(DynamicCast<SolidParticles>(this, body.base_particles_)),
		  solid_body_density_(DynamicCast<Solid>(this, base_particles_->base_material_)->ReferenceDensity())
	{
		setMassProperties();
	}
	//=================================================================================================//
}
