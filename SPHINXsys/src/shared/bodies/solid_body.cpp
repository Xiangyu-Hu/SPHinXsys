/**
 * @file    solid_body.cpp
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "solid_body.h"
#include "mesh_cell_linked_list.h"
#include "elastic_solid.h"
#include "solid_particles.h"
#include "sph_system.h"
#include "base_kernel.h"

namespace SPH {
	//===============================================================//
	SolidBody::SolidBody(SPHSystem &system, string body_name, Solid &solid_material, 
		int refinement_level, ParticlesGeneratorOps op)
		: RealBody(system, body_name, solid_material, refinement_level, 1.0, op)
	{
	
	}
	//===============================================================//
	void SolidBody::BuildInnerConfiguration()
	{
		mesh_cell_linked_list_->BuildReferenceInnerConfiguration(*this);
	}
	//===============================================================//
	void SolidBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->BuildReferenceContactConfiguration(*this);
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===============================================================//
	SolidBodyPart::SolidBodyPart(SolidBody *solid_body, string soild_body_part_name)
		: LagrangianBodyPart(solid_body, soild_body_part_name),
		solid_body_(solid_body), soild_body_part_region_(soild_body_part_name)
	{

	}
	//===============================================================//
	void SolidBodyPart::TagBodyPartParticles()
	{
		for (size_t i = 0; i < solid_body_->number_of_particles_; ++i)
		{
			BaseParticleData &base_particle_data_i
				= solid_body_->base_particles_->base_particle_data_[i];

			if (soild_body_part_region_.contain(base_particle_data_i.pos_n_))
			{
				body_part_particles_.push_back(i);
			}
		}
	}
	//===============================================================//
	SolidBodyPartForSimbody
		::SolidBodyPartForSimbody(SolidBody *solid_body,
			string soild_body_part_name, Real solid_body_density)
		: SolidBodyPart(solid_body, soild_body_part_name), 
		solid_body_density_(solid_body_density)
	{

	}
}