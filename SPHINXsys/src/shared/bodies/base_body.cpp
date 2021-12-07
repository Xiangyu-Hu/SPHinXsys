/**
 * @file 	base_body.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	Chi ZHang and Xiangyu Hu
 */
#include "base_body.h"

#include "sph_system.h"
#include "base_particles.h"
#include "body_relation.h"
#include "cell_linked_list.h"

namespace SPH
{
	//=================================================================================================//
	SPHBody::SPHBody(SPHSystem &sph_system, std::string body_name,
					 ParticleAdaptation *particle_adaptation, ParticleGenerator *particle_generator)
		: sph_system_(sph_system), body_name_(body_name), newly_updated_(true),
		  body_domain_bounds_(0, 0), is_domain_bounds_determined_(false),
		  particle_adaptation_(particle_adaptation), particle_generator_(particle_generator),
		  body_shape_(nullptr), generative_structure_(nullptr)
	{
		sph_system_.addABody(this);
		particle_adaptation_->initialize(this);
	}
	//=================================================================================================//
	BoundingBox SPHBody::getSPHSystemBounds()
	{
		return sph_system_.system_domain_bounds_;
	}
	//=================================================================================================//
	std::string SPHBody::getBodyName()
	{
		return body_name_;
	}
	//=================================================================================================//
	SPHSystem &SPHBody::getSPHSystem()
	{
		return sph_system_;
	}
	//=================================================================================================//
	void SPHBody::useParticleGeneratorReload()
	{
		particle_generator_->~ParticleGenerator();
		particle_generator_ = new ParticleGeneratorReload(sph_system_.in_output_, body_name_);
	}
	//=================================================================================================//
	void SPHBody::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
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
	void SPHBody::setBodyDomainBounds(BoundingBox body_domain_bounds) 
	{
		body_domain_bounds_ = body_domain_bounds;
		is_domain_bounds_determined_ = true;
	};
	//=================================================================================================//
	BoundingBox SPHBody::getBodyDomainBounds()
	{
		if (!is_domain_bounds_determined_)
		{
			body_domain_bounds_ = body_shape_->findBounds();
			is_domain_bounds_determined_ = true;
		}
		return body_domain_bounds_;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToVtuFile(std::ostream &output_file)
	{
		base_particles_->writeParticlesToVtuFile(output_file);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeSurfaceParticlesToVtuFile(std::ostream &output_file, ShapeSurface& surface_particles)
	{
		base_particles_->writeSurfaceParticlesToVtuFile(output_file, surface_particles);
		newly_updated_ = false;
	}
	//=================================================================================================//
	void SPHBody::writeParticlesToPltFile(std::ofstream &output_file)
	{
		if (newly_updated_)
			base_particles_->writeParticlesToPltFile(output_file);
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
	RealBody::RealBody(SPHSystem &sph_system, std::string body_name,
					   ParticleAdaptation *particle_adaptation, ParticleGenerator *particle_generator)
		: SPHBody(sph_system, body_name, particle_adaptation, particle_generator),
		  particle_sorting_(this)
	{
		sph_system.addARealBody(this);
		cell_linked_list_ = particle_adaptation_->createCellLinkedList();
		size_t number_of_split_cell_lists = powerN(3, Vecd(0).size());
		split_cell_lists_.resize(number_of_split_cell_lists);
	}
	//=================================================================================================//
	void RealBody::assignBaseParticles(BaseParticles *base_particles)
	{
		SPHBody::assignBaseParticles(base_particles);
		particle_sorting_.assignBaseParticles(base_particles);
		cell_linked_list_->assignBaseParticles(base_particles);
	}
	//=================================================================================================//
	void RealBody::sortParticleWithCellLinkedList()
	{
		StdLargeVec<size_t> &sequence = base_particles_->sequence_;
		size_t size = base_particles_->total_real_particles_;
		cell_linked_list_->computingSequence(sequence);
		particle_sorting_.sortingParticleData(sequence.data(), size);
	}
	//=================================================================================================//
	void RealBody::updateCellLinkedList()
	{
		cell_linked_list_->UpdateCellLists();
	}
	FictitiousBody::
		FictitiousBody(SPHSystem &system, std::string body_name,
					   ParticleAdaptation *particle_adaptation, ParticleGenerator *particle_generator)
		: SPHBody(system, body_name, particle_adaptation, particle_generator)
	{
		system.addAFictitiousBody(this);
	}
	//=================================================================================================//
	BodyPartByShape::BodyPartByShape(SPHBody *body, std::string body_part_name)
		: BodyPart(body, body_part_name), body_part_shape_(nullptr) {}
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
		BaseParticles *base_particles = body_->base_particles_;
		for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
		{
			if (body_part_shape_->checkContain(base_particles->pos_n_[i]))
				tagAParticle(i);
		}
	}
	//=================================================================================================//
	ShapeSurface::ShapeSurface(SPHBody *body)
		: BodyPartByParticle(body, "Surface"),
		  particle_spacing_min_(body->particle_adaptation_->MinimumSpacing())
	{
		tagBodyPart();
	}
	//=================================================================================================//
	void ShapeSurface::tagBodyPart()
	{
		BaseParticles *base_particles = body_->base_particles_;
		for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
		{
			Real phi = body_->body_shape_->findSignedDistance(base_particles->pos_n_[i]);
			if (fabs(phi) < particle_spacing_min_)
				tagAParticle(i);
		}
		//std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
	}
	//=================================================================================================//
	ShapeSurfaceLayer::ShapeSurfaceLayer(SPHBody *body, Real layer_thickness)
		: BodyPartByParticle(body, "InnerLayers"),
		  thickness_threshold_(body->particle_adaptation_->ReferenceSpacing() * layer_thickness)
	{
		tagBodyPart();
	}
	//=================================================================================================//
	void ShapeSurfaceLayer::tagBodyPart()
	{
		BaseParticles *base_particles = body_->base_particles_;
		for (size_t i = 0; i < base_particles->total_real_particles_; ++i)
		{
			Vecd position_i = base_particles->pos_n_[i];
			Real distance = fabs(body_->body_shape_->findSignedDistance(position_i));
			if (distance < thickness_threshold_)
				tagAParticle(i);
		}
		std::cout << "Number of inner layers particles : " << body_part_particles_.size() << std::endl;
	}
	//=================================================================================================//
	BodyPartByCell::BodyPartByCell(RealBody *real_body, std::string body_part_name)
		: BodyPartByShape(real_body, body_part_name), real_body_(real_body),
		  checkIncluded_(std::bind(&BodyPartByCell::checkIncluded, this, _1, _2)) {}
	//=================================================================================================//
	bool BodyPartByCell::checkIncluded(Vecd cell_position, Real threshold)
	{
		return body_part_shape_->checkNotFar(cell_position, threshold);
	}
	//=================================================================================================//
	void BodyPartByCell::tagBodyPart()
	{
		real_body_->cell_linked_list_->tagBodyPartByCell(body_part_cells_, checkIncluded_);
	}
	//=================================================================================================//
	NearShapeSurface::
		NearShapeSurface(RealBody *real_body, ComplexShape *complex_shape, std::string body_part_name)
		: BodyPartByCell(real_body, body_part_name)
	{
		level_set_complex_shape_ = new LevelSetComplexShape(real_body, *complex_shape, true);
		body_part_shape_ = level_set_complex_shape_;
		tagBodyPart();
	}
	//=================================================================================================//
	NearShapeSurface::NearShapeSurface(RealBody *real_body)
		: BodyPartByCell(real_body, "NearShapeSurface")
	{
		body_part_shape_ = real_body->body_shape_;
		level_set_complex_shape_ = dynamic_cast<LevelSetComplexShape *>(body_part_shape_);
		if (level_set_complex_shape_ == nullptr)
		{
			std::cout << "\n FAILURE: LevelSetComplexShape is undefined!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		tagBodyPart();
	}
	//=================================================================================================//
	LevelSetComplexShape *NearShapeSurface::getLevelSetComplexShape()
	{
		return level_set_complex_shape_;
	}
	//=================================================================================================//
	bool NearShapeSurface::checkIncluded(Vecd cell_position, Real threshold)
	{
		return body_part_shape_->checkNearSurface(cell_position, threshold);
	}
	//=================================================================================================//
	TerminateBranches::TerminateBranches(SPHBody *body) 
		: BodyPartByParticle(body, "Leaves"),
		  tree_(dynamic_cast<GenerativeTree *>(body->generative_structure_))
	{
		tagBodyPart();
	}
	//=================================================================================================//
	void TerminateBranches::tagBodyPart()
	{
		for (const auto	*branch : tree_->branches_)
		{
			if (branch->is_terminated_)
			{
				size_t particle_id = branch->inner_particles_.back();
				tagAParticle(particle_id);
			}
		}
	}
	//=================================================================================================//
}
