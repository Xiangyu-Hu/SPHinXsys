/**
 * @file 	base_body.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	hi ZHang and Xiangyu Hu
 * @version	0.1
 * 			0.2.0
 * 			Cell splitting algorithm are added. 
 * 			Chi Zhang
 */
#include "base_body.h"
#include "sph_system.h"
#include "in_output.h"
#include "base_particles.h"
#include "all_kernels.h"
#include "mesh_cell_linked_list.h"
//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	SPHBody::SPHBody(SPHSystem &sph_system, string body_name,
		int refinement_level, Real smoothing_length_ratio, ParticlesGeneratorOps op)
	: sph_system_(sph_system), body_shape_(body_name), body_name_(body_name), newly_updated_(true),
		refinement_level_(refinement_level), particle_generator_op_(op),
		body_lower_bound_(0), body_upper_bound_(0), prescribed_body_bounds_(false),
		levelset_mesh_(NULL)
	{	
		sph_system_.AddBody(this);

		particle_spacing_ 	= RefinementLevelToParticleSpacing();
		smoothing_length_ = particle_spacing_ * smoothing_length_ratio;
		kernel_ 			= GenerateAKernel(smoothing_length_);
		base_mesh_cell_linked_list_
							= new MeshCellLinkedList(this, sph_system.lower_bound_,
									sph_system_.upper_bound_, kernel_->GetCutOffRadius());
		size_t number_of_split_cell_lists = powern(3, Vecd(0).size());
		/** I will use concurrent vector here later after tests. */
		split_cell_lists_.resize(number_of_split_cell_lists);
	}
	//=================================================================================================//
	void  SPHBody::getSPHSystemBound(Vecd& system_lower_bound, Vecd& system_uppwer_bound) 
	{
		system_lower_bound = sph_system_.lower_bound_;
		system_uppwer_bound = sph_system_.upper_bound_;
	}
	//=================================================================================================//
	Real SPHBody::RefinementLevelToParticleSpacing()
	{
		return sph_system_.particle_spacing_ref_	
			/ powern(2.0, refinement_level_);
	}
	//=================================================================================================//
	Kernel* SPHBody::GenerateAKernel(Real smoothing_length)
	{
		return new KernelWendlandC2(smoothing_length);
	}
	//=================================================================================================//
	void SPHBody::ReplaceKernelFunction(Kernel* kernel)
	{
		delete kernel_;
		kernel_ = kernel;
		base_mesh_cell_linked_list_->reassignKernel(kernel);
	}
	//=================================================================================================//
	string SPHBody::GetBodyName()
	{
		return body_name_;
	}
	//=================================================================================================//
	void SPHBody::AllocateConfigurationMemoriesForBodyBuffer()
	{
		for (size_t i = 0; i < body_relations_.size(); i++)
		{
			body_relations_[i]->updateConfigurationMemories();
		}
	}
	//=================================================================================================//
	void SPHBody::addLevelsetMesh(Real mesh_size_ratio)
	{
		Vecd body_lower_bound, body_upper_bound;
		findBodyShapeBounds(body_lower_bound, body_upper_bound);
		/** levelset mesh. */
		Real mesh_spacing = mesh_size_ratio * particle_spacing_;

		levelset_mesh_ = new LevelSet(
			this, 				/**< body pointer. */
			body_lower_bound, 	/**< Lower bound. */
			body_upper_bound,	/**< Upper bound. */
			mesh_spacing, 		/**< Mesh spacing. */
			4					/**< Buffer size. */
		);
	}
	//=================================================================================================//

	bool SPHBody::checkBodyShapeContain(Vecd pnt)
	{
		return body_shape_.checkContain(pnt);
	}
	//=================================================================================================//
	Vecd SPHBody::ClosestPointOnBodySurface(Vecd input_pnt)
	{
		return body_shape_.findClosestPoint(input_pnt);
	}
	//=================================================================================================//
	void SPHBody::findBodyShapeBounds(Vecd& lower_bound, Vecd& upper_bound)
	{
		if(!prescribed_body_bounds_) 
		{
			body_shape_.findBounds(lower_bound, upper_bound);
		}
		else 
		{
			lower_bound = body_lower_bound_;
			upper_bound = body_upper_bound_;
		}
	}
	//=================================================================================================//
	void SPHBody::WriteParticlesToVtuFile(ofstream &output_file)
	{
		base_particles_->WriteParticlesToVtuFile(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::WriteParticlesToPltFile(ofstream &output_file)
	{
		if (newly_updated_) base_particles_->WriteParticlesToPltFile(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::WriteParticlesToXmlForRestart(std::string &filefullpath)
	{
		base_particles_->WriteParticlesToXmlForRestart(filefullpath);
	}
	//=================================================================================================//
	void SPHBody::ReadParticlesFromXmlForRestart(std::string &filefullpath)
	{
		base_particles_->ReadParticleFromXmlForRestart(filefullpath);
	}
	//=================================================================================================//
	void SPHBody::WriteToXmlForReloadParticle(std::string &filefullpath)
	{
		base_particles_->WriteToXmlForReloadParticle(filefullpath);
	}
	//=================================================================================================//
	void SPHBody::ReadFromXmlForReloadParticle(std::string &filefullpath)
	{
		base_particles_->ReadFromXmlForReloadParticle(filefullpath);
	}
	//=================================================================================================//
	SPHBody* SPHBody::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	RealBody::RealBody(SPHSystem &sph_system, string body_name,
		int refinement_level, Real smoothing_length_ratio, ParticlesGeneratorOps op)
	: SPHBody(sph_system, body_name, refinement_level, smoothing_length_ratio, op)
	{
		sph_system.AddRealBody(this);

		base_mesh_cell_linked_list_->allocateMeshDataMatrix();
	}
	//=================================================================================================//
	void RealBody::AllocateMemoryCellLinkedList()
	{
		base_mesh_cell_linked_list_->allocateMeshDataMatrix();
	}
	//=================================================================================================//
	void RealBody::UpdateCellLinkedList()
	{
		base_mesh_cell_linked_list_->UpdateCellLists();
	}
	//=================================================================================================//
	RealBody* RealBody::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	FictitiousBody::FictitiousBody(SPHSystem &system, string body_name, int refinement_level,
		Real smoothing_length_ratio, ParticlesGeneratorOps op)
	: SPHBody(system, body_name, refinement_level, smoothing_length_ratio, op)
	{
		system.AddFictitiousBody(this);
	}
	//=================================================================================================//
	void FictitiousBody::AllocateMemoryCellLinkedList()
	{
		/** do nothing here. */;
	}
	//=================================================================================================//
	void FictitiousBody::UpdateCellLinkedList()
	{
		/** do nothing here. */;
	}
	//=================================================================================================//
	FictitiousBody* FictitiousBody::PointToThisObject()
	{
		return this;
	}
	//=================================================================================================//
	void BodyPartByParticle::tagAParticle(size_t particle_index)
	{
		BaseParticleData& base_particle_data_i
			= body_->base_particles_->base_particle_data_[particle_index];
		body_part_particles_.push_back(particle_index);
		base_particle_data_i.is_sortable_ = false;
	}
	//=================================================================================================//
	void BodyPartByParticle::TagBodyPart()
	{
		BaseParticles* base_particles = body_->base_particles_;
		for (size_t i = 0; i < body_->number_of_particles_; ++i)
		{
			BaseParticleData &base_particle_data_i
				= base_particles->base_particle_data_[i];
			if (body_part_shape_.checkContain(base_particle_data_i.pos_n_)) tagAParticle(i);
		}
	}
	//=================================================================================================//
	BodySurface::BodySurface(SPHBody* body)
		: BodyPartByParticle(body, "Surface")
	{
		if (body->levelset_mesh_ == NULL)
		{
			std::cout << "\n BodySurface::BodySurface: Background mesh is required. Exit the program! \n";
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(0);
		}
		TagBodyPart();
	}	
	//=================================================================================================//
	void BodySurface::TagBodyPart()
	{
		for (size_t i = 0; i < body_->number_of_particles_; ++i)
		{
			BaseParticleData& base_particle_data_i
				= body_->base_particles_->base_particle_data_[i];

			Real phii = body_->levelset_mesh_->probeLevelSet(base_particle_data_i.pos_n_);
			//this is important, as outer particles is neglect, is should be the particle spacing
			if (phii < body_->particle_spacing_) tagAParticle(i);
		}
		std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
	}
	//=================================================================================================//
	BodySurfaceLayer::BodySurfaceLayer(SPHBody* body, Real layer_thickness)
		: BodyPartByParticle(body, "InnerLayers"), layer_thickness_(layer_thickness)
	{
		TagBodyPart();
	}
	//=================================================================================================//
	void BodySurfaceLayer::TagBodyPart()
	{
		for (size_t i = 0; i < body_->number_of_particles_; ++i)
		{
			BaseParticleData& base_particle_data_i
				= body_->base_particles_->base_particle_data_[i];
			Real distance 
				= (body_->ClosestPointOnBodySurface(base_particle_data_i.pos_n_) - base_particle_data_i.pos_n_).norm();
			if (distance < body_->particle_spacing_* layer_thickness_) tagAParticle(i);
		}
		std::cout << "Number of inner layers particles : " << body_part_particles_.size() << std::endl;
	}
	//=================================================================================================//
	NearBodySurface::NearBodySurface(SPHBody* body)
		: BodyPartByCell(body, "NearBodySurface")
	{
		if (body->levelset_mesh_ == NULL)
		{
			std::cout << "\n NearBodySurface::NearBodySurface: Levelset mesh is required. Exit the program! \n";
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(0);
		}
		TagBodyPart();
	}
	//=================================================================================================//
	SolidBodyPartForSimbody
		::SolidBodyPartForSimbody(SPHBody* solid_body, string solid_body_part_name)
		: BodyPartByParticle(solid_body, solid_body_part_name)
	{
		Solid* solid = dynamic_cast<Solid*>(body_->base_particles_->base_material_);
		solid_body_density_ = solid->getReferenceDensity();
	}
	//=================================================================================================//
}