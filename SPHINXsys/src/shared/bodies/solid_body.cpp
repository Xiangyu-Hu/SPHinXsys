/**
 * @file    solid_body.cpp
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "solid_body.h"
#include "sph_system.h"
#include "mesh_cell_linked_list.h"

namespace SPH {
	//===============================================================//
	SolidBody::SolidBody(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: RealBody(system, body_name, refinement_level, 1.3, op)
	{
	
	}
	//===============================================================//
	void SolidBody::BuildInnerConfiguration()
	{
		mesh_cell_linked_list_->UpdateInnerConfiguration(*this, reference_inner_configuration_);
	}
	//===============================================================//
	void SolidBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===============================================================//
	SolidBodyPartForSimbody
		::SolidBodyPartForSimbody(SolidBody *solid_body,
			string soild_body_part_name)
		: BodyPartByParticle(solid_body, soild_body_part_name)
	{
		Solid* solid = dynamic_cast<Solid*>(body_->base_particles_->base_material_);
		solid_body_density_ = solid->getReferenceDensity();
	}
}