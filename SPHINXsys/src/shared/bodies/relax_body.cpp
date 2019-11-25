/**
 * @file 	relax_body.cpp
 * @brief 	This is the class for bodies used for particle relaxation scheme.
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "relax_body.h"
#include "mesh_cell_linked_list.h"
#include "sph_system.h"
#include "base_kernel.h"
#include "relax_body_particles.h"
#include "base_mesh.h"

namespace SPH {
	//===========================================================//
	RelaxBody::RelaxBody(SPHSystem &sph_system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: RealBody(sph_system, body_name, refinement_level, 1.0, op)
	{

	}
	//===========================================================//
	void RelaxBody::InitializeBackgroundMesh()
	{
		Vecd body_lower_bound, body_upper_bound;
		BodyBounds(body_lower_bound, body_upper_bound);
		//background mesh has much higher resolution
		mesh_background_
			= new MeshBackground(body_lower_bound,
				body_upper_bound, 0.5*particle_spacing_, 4);
		mesh_background_->AllocateMeshDataMatrix();
		mesh_background_->InitializeLevelSetData(*this);
		mesh_background_->ComputeCurvatureFromLevelSet(*this);
	}
	//===============================================================//
	void RelaxBody::BuildInnerConfiguration()
	{
		mesh_cell_linked_list_->UpdateInnerConfiguration(*this);
	}
	//===============================================================//
	void RelaxBody::BuildContactConfiguration()
	{
		mesh_cell_linked_list_->UpdateContactConfiguration(*this);
	}
	//===========================================================//
	RelaxBodySurface::RelaxBodySurface(RelaxBody *relax_body)
		: BodyPartByParticle(relax_body, "Surface"), relax_body_(relax_body)
	{
		TagBodyPartParticles();
	}
	//===========================================================//
	void RelaxBodySurface::TagBodyPartParticles()
	{
		for (size_t i = 0; i < relax_body_->number_of_particles_; ++i)
		{
			BaseParticleData &base_particle_data_i
				= relax_body_->base_particles_->base_particle_data_[i];

			Real phii = relax_body_->mesh_background_->ProbeLevelSet(base_particle_data_i.pos_n_);
			//this is important, as outer particles is neglect, is shoul be the particle spacing
			if (phii < relax_body_->particle_spacing_)
			{
				body_part_particles_.push_back(i);
			}
		}
		std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
	}
	//===========================================================//
}