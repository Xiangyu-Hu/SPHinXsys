/**
 * @file    solid_body.cpp
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "solid_body.h"

#include "sph_system.h"
#include "base_material.h"
#include "solid_particles.h"

namespace SPH
{
	//=================================================================================================//
	SolidBody::SolidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: RealBody(system, shape_ptr)
	{
		sph_system_.solid_bodies_.push_back(this);
		defineAdaptation<SPHAdaptation>(1.15);
	}
	//=================================================================================================//
	SolidBodyPartForSimbody::
		SolidBodyPartForSimbody(SPHBody &body,SharedPtr<Shape> shape_ptr)
		: BodyRegionByParticle(body, shape_ptr),
		  solid_particles_(DynamicCast<SolidParticles>(this, &body.getBaseParticles())),
		  solid_body_density_(DynamicCast<Solid>(this, body.base_material_)->ReferenceDensity())
	{
		setMassProperties();
	}
	//=================================================================================================//
}
