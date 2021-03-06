/**
 * @file 	base_body.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	Chi ZHang and Xiangyu Hu
 */
#include "base_body.h"

#include "sph_system.h"
#include "in_output.h"
#include "base_particles.h"
#include "mesh_cell_linked_list.h"
#include "geometry_level_set.h"

//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	SPHBody::SPHBody(SPHSystem &sph_system, string body_name,
		ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator) :
		sph_system_(sph_system), body_name_(body_name), newly_updated_(true),
		body_domain_bounds_(0, 0), prescribed_body_bounds_(false),
		particle_adaptation_(particle_adaptation), particle_generator_(particle_generator),
		body_shape_(NULL)
	{	
		sph_system_.addABody(this);
		particle_adaptation_->initialize(this);
	}
	//=================================================================================================//
	BoundingBox  SPHBody::getSPHSystemBounds()
	{
		return sph_system_.system_domain_bounds_;
	}
	//=================================================================================================//
	string SPHBody::getBodyName()
	{
		return body_name_;
	}
	//=================================================================================================//
	SPHSystem& SPHBody::getSPHSystem()
	{
		return sph_system_;
	}
	//=================================================================================================//
	void SPHBody::assignBaseParticles(BaseParticles* base_particles)
	{
		base_particles_ = base_particles;
	}
	//=================================================================================================//
	Real SPHBody::computeReferenceNumberDensity(Vec2d zero)
	{
		Real sigma(0);
		Kernel* kernel = particle_adaptation_->getKernel();
		Real cutoff_radius = kernel->CutOffRadius();
		Real particle_spacing = particle_adaptation_->ReferenceSpacing();
		int search_range = int(cutoff_radius / particle_spacing) + 1;
		for (int j = -search_range; j <= search_range; ++j)
			for (int i = -search_range; i <= search_range; ++i)
			{
				Vec2d particle_location(Real(i) * particle_spacing, Real(j) * particle_spacing);
				Real distance = particle_location.norm();
				if (distance < cutoff_radius) sigma += kernel->W(distance, particle_location);
			}
		return sigma;
	}
	//=================================================================================================//
	Real SPHBody::computeReferenceNumberDensity(Vec3d zero)
	{
		Real sigma(0);
		Kernel* kernel = particle_adaptation_->getKernel();
		Real cutoff_radius = kernel->CutOffRadius();
		Real particle_spacing = particle_adaptation_->ReferenceSpacing();
		int search_range = int(cutoff_radius / particle_spacing) + 1;
		for (int k = -search_range; k <= search_range; ++k)
			for (int j = -search_range; j <= search_range; ++j)
				for (int i = -search_range; i <= search_range; ++i)
				{
					Vec3d particle_location(Real(i) * particle_spacing,
						Real(j) * particle_spacing, Real(k) * particle_spacing);
					Real distance = particle_location.norm();
					if (distance < cutoff_radius) sigma += kernel->W(distance, particle_location);
				}
		return sigma;
	}
	//=================================================================================================//
	void SPHBody::allocateConfigurationMemoriesForBufferParticles()
	{
		for (size_t i = 0; i < body_relations_.size(); i++)
		{
			body_relations_[i]->updateConfigurationMemories();
		}
	}
	//=================================================================================================//
	BoundingBox SPHBody::findBodyDomainBounds()
	{
		if(!prescribed_body_bounds_) 
		{
			return body_shape_->findBounds();
		}
		else 
		{
			return body_domain_bounds_;
		}
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToVtuFile(ofstream &output_file)
	{
		base_particles_->writeParticlesToVtuFile(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToPltFile(ofstream &output_file)
	{
		if (newly_updated_) base_particles_->writeParticlesToPltFile(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		base_particles_->writeParticlesToXmlForRestart(filefullpath);
	}
	//=================================================================================================//
	void SPHBody::readParticlesFromXmlForRestart(std::string &filefullpath)
	{
		base_particles_->readParticleFromXmlForRestart(filefullpath);
	}
	//=================================================================================================//
	void SPHBody::writeToXmlForReloadParticle(std::string &filefullpath)
	{
		base_particles_->writeToXmlForReloadParticle(filefullpath);
	}
	//=================================================================================================//
	void SPHBody::readFromXmlForReloadParticle(std::string &filefullpath)
	{
		base_particles_->readFromXmlForReloadParticle(filefullpath);
	}
	//=================================================================================================//
	RealBody::RealBody(SPHSystem &sph_system, string body_name,
		ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator) : 
		SPHBody(sph_system, body_name, particle_adaptation, particle_generator),
		particle_sorting_(this)
	{
		sph_system.addARealBody(this);
		mesh_cell_linked_list_ = particle_adaptation_->createMeshCellLinkedList();
		size_t number_of_split_cell_lists = powern(3, Vecd(0).size());
		split_cell_lists_.resize(number_of_split_cell_lists);
	}
	//=================================================================================================//
	void RealBody::assignBaseParticles(BaseParticles* base_particles)
	{
		SPHBody::assignBaseParticles(base_particles);
		particle_sorting_.assignBaseParticles(base_particles);
		mesh_cell_linked_list_->assignBaseParticles(base_particles);
	}
	//=================================================================================================//
	void RealBody::sortParticleWithMeshCellLinkedList()
	{
		StdLargeVec<size_t>& sequence = base_particles_->sequence_;
		size_t size = base_particles_->total_real_particles_;
		mesh_cell_linked_list_->computingSequence(sequence);
		particle_sorting_.sortingParticleData(sequence.data(), size);
	}
	//=================================================================================================//
	void RealBody::updateCellLinkedList()
	{
		mesh_cell_linked_list_->UpdateCellLists();
	}
	FictitiousBody::FictitiousBody(SPHSystem &system, string body_name, 
		ParticleAdaptation* particle_adaptation, ParticleGenerator* particle_generator)
	: SPHBody(system, body_name, particle_adaptation, particle_generator)
	{
		system.addAFictitiousBody(this);
	}
	//=================================================================================================//
	BodyPartByShape::BodyPartByShape(SPHBody* body, string body_part_name) : 
	BodyPart(body, body_part_name), body_part_shape_(NULL) {}	
	//=================================================================================================//
	BoundingBox BodyPartByShape::BodyPartBounds() 
	{ 
		return body_part_shape_->findBounds(); 
	}
	//=================================================================================================//
	void BodyPartByParticle::tagAParticle(size_t particle_index)
	{
		body_part_particles_.push_back(particle_index);
	}
	//=================================================================================================//
	void BodyPartByParticle::tagBodyPart()
	{
		BaseParticles* base_particles = body_->base_particles_;
		for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
		{
			if (body_part_shape_->checkContain(base_particles->pos_n_[i])) tagAParticle(i);
		}
	}
	//=================================================================================================//
	ShapeSurface::ShapeSurface(SPHBody* body)
		: BodyPartByParticle(body, "Surface"),
		particle_spacing_min_(body->particle_adaptation_->MinimumSpacing())
	{
		tagBodyPart();
	}	
	//=================================================================================================//
	void ShapeSurface::tagBodyPart()
	{
		BaseParticles* base_particles = body_->base_particles_;
		for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
		{
			Real phi = body_->body_shape_->findSignedDistance(base_particles->pos_n_[i]);
			if (fabs(phi) < particle_spacing_min_) tagAParticle(i);
		}
		std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
	}
	//=================================================================================================//
	ShapeSurfaceLayer::ShapeSurfaceLayer(SPHBody* body, Real layer_thickness)
		: BodyPartByParticle(body, "InnerLayers"), 
		thickness_threshold_(body->particle_adaptation_->ReferenceSpacing() * layer_thickness)
	{
		tagBodyPart();
	}
	//=================================================================================================//
	void ShapeSurfaceLayer::tagBodyPart()
	{
		BaseParticles* base_particles = body_->base_particles_;
		for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
		{
			Vecd position_i = base_particles->pos_n_[i];
			Real distance = fabs(body_->body_shape_->findSignedDistance(position_i));
			if (distance < thickness_threshold_) tagAParticle(i);
		}
		std::cout << "Number of inner layers particles : " << body_part_particles_.size() << std::endl;
	}
	//=================================================================================================//
	 BodyPartByCell::BodyPartByCell(RealBody *real_body, string body_part_name)	: 
	 	BodyPartByShape(real_body, body_part_name), real_body_(real_body),
		checkIncluded_(std::bind(&BodyPartByCell::checkIncluded, this, _1, _2)) {}
	//=================================================================================================//
	bool BodyPartByCell::checkIncluded(Vecd cell_position, Real threshold)
	{
		return body_part_shape_->checkNotFar(cell_position, threshold);
	}
	//=================================================================================================//
	void BodyPartByCell::tagBodyPart()
	{
		real_body_->mesh_cell_linked_list_->tagBodyPartByCell(body_part_cells_, checkIncluded_);
	}
	//=================================================================================================//
	NearShapeSurface::NearShapeSurface(RealBody* real_body, ComplexShape* complex_shape, string body_part_name) :
		BodyPartByCell(real_body, body_part_name)
	{
		level_set_complex_shape_ = new LevelSetComplexShape(real_body, *complex_shape, true);
		body_part_shape_ = level_set_complex_shape_;
		tagBodyPart();
	}
	//=================================================================================================//
	NearShapeSurface::NearShapeSurface(RealBody* real_body) :
		BodyPartByCell(real_body, "NearShapeSurface")
	{
		body_part_shape_ = real_body->body_shape_;
		level_set_complex_shape_ = dynamic_cast<LevelSetComplexShape*>(body_part_shape_);
		if (level_set_complex_shape_ == NULL)
		{
			std::cout << "\n FAILURE: LevelSetComplexShape is undefined!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		tagBodyPart();
	}
	//=================================================================================================//
	LevelSetComplexShape* NearShapeSurface::getLevelSetComplexShape() 
	{
		return level_set_complex_shape_;
	}	
	//=================================================================================================//
	bool NearShapeSurface::checkIncluded(Vecd cell_position, Real threshold)
	{
		return body_part_shape_->checkNearSurface(cell_position, threshold);
	}	
	//=================================================================================================//
}
